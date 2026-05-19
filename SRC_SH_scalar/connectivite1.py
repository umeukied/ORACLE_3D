import re
from collections import defaultdict

VTK_FILE = "gortlercosin.vtk"


def read_legacy_vtk_unstructured_grid(filename):
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # -----------------------------
    # POINTS
    # -----------------------------
    points = []
    i = 0
    while i < len(lines):
        if lines[i].startswith("POINTS"):
            parts = lines[i].split()
            npoints = int(parts[1])
            i += 1

            coords = []
            while len(coords) < 3 * npoints:
                coords.extend(lines[i].split())
                i += 1

            coords = list(map(float, coords))
            for k in range(npoints):
                x = coords[3 * k]
                y = coords[3 * k + 1]
                z = coords[3 * k + 2]
                points.append((x, y, z))
            break
        i += 1

    # -----------------------------
    # CELLS
    # -----------------------------
    cells = []
    while i < len(lines):
        if lines[i].startswith("CELLS"):
            parts = lines[i].split()
            ncells = int(parts[1])
            i += 1

            for _ in range(ncells):
                vals = list(map(int, lines[i].split()))
                nverts = vals[0]
                conn = vals[1:]
                cells.append(conn)
                i += 1
            break
        i += 1

    # -----------------------------
    # CELL_TYPES
    # -----------------------------
    cell_types = []
    while i < len(lines):
        if lines[i].startswith("CELL_TYPES"):
            parts = lines[i].split()
            ncelltypes = int(parts[1])
            i += 1

            while len(cell_types) < ncelltypes:
                cell_types.extend(map(int, lines[i].split()))
                i += 1
            break
        i += 1

    return points, cells, cell_types


def build_triangle_neighbors(cells, cell_types):
    """
    Ne garde que les triangles VTK type 5.
    Retourne:
      triangles: liste des triangles [(p0,p1,p2), ...]
      tri_to_original: indice triangle -> indice cellule VTK originale
      neighbors: pour chaque triangle, voisins par face
    """
    triangles = []
    tri_to_original = []

    for icell, (conn, ctype) in enumerate(zip(cells, cell_types)):
        if ctype == 5 and len(conn) == 3:  # triangle
            triangles.append(tuple(conn))
            tri_to_original.append(icell)

    # face -> liste des triangles qui portent cette face
    face_to_tris = defaultdict(list)

    for itri, (a, b, c) in enumerate(triangles):
        faces = [
            tuple(sorted((a, b))),
            tuple(sorted((b, c))),
            tuple(sorted((c, a))),
        ]
        for iface, face in enumerate(faces):
            face_to_tris[face].append((itri, iface))

    neighbors = []
    for itri, (a, b, c) in enumerate(triangles):
        tri_faces = [
            (a, b),
            (b, c),
            (c, a),
        ]

        tri_neighbors = []
        for face in tri_faces:
            fkey = tuple(sorted(face))
            attached = face_to_tris[fkey]

            # S'il n'y a qu'un triangle sur cette face => bord
            if len(attached) == 1:
                tri_neighbors.append(-1)
            else:
                # Chercher l'autre triangle
                other = [t for t, _ in attached if t != itri]
                tri_neighbors.append(other[0] if other else -1)

        neighbors.append(tri_neighbors)

    return triangles, tri_to_original, neighbors


def export_triangle_info(points, triangles, neighbors, out_file="mesh_info.txt"):
    with open(out_file, "w") as f:
        for ielem, tri in enumerate(triangles):
            a, b, c = tri
            xa, ya, za = points[a]
            xb, yb, zb = points[b]
            xc, yc, zc = points[c]

            f.write(f"Element {ielem}\n")
            f.write(f"  Type           : triangle\n")
            f.write(f"  Number of faces: 3\n")
            f.write(f"  Vertex IDs     : {a}, {b}, {c}\n")
            f.write(f"  Vertex coords  :\n")
            f.write(f"    P{a} = ({xa}, {ya}, {za})\n")
            f.write(f"    P{b} = ({xb}, {yb}, {zb})\n")
            f.write(f"    P{c} = ({xc}, {yc}, {zc})\n")

            faces = [(a, b), (b, c), (c, a)]
            for iface, (face, neigh) in enumerate(zip(faces, neighbors[ielem]), start=1):
                f.write(f"  Face {iface}       : {face}, neighbor = {neigh}\n")

            f.write("\n")


def export_csv(points, triangles, neighbors, out_file="mesh_info.csv"):
    with open(out_file, "w") as f:
        f.write(
            "elem_id,p1,p2,p3,"
            "x1,y1,z1,x2,y2,z2,x3,y3,z3,"
            "nfaces,neighbor_face1,neighbor_face2,neighbor_face3\n"
        )

        for ielem, tri in enumerate(triangles):
            a, b, c = tri
            xa, ya, za = points[a]
            xb, yb, zb = points[b]
            xc, yc, zc = points[c]
            n1, n2, n3 = neighbors[ielem]

            f.write(
                f"{ielem},{a},{b},{c},"
                f"{xa},{ya},{za},{xb},{yb},{zb},{xc},{yc},{zc},"
                f"3,{n1},{n2},{n3}\n"
            )


def main():
    points, cells, cell_types = read_legacy_vtk_unstructured_grid(VTK_FILE)

    triangles, tri_to_original, neighbors = build_triangle_neighbors(cells, cell_types)

    print(f"Nombre de points      : {len(points)}")
    print(f"Nombre de cellules    : {len(cells)}")
    print(f"Nombre de triangles   : {len(triangles)}")

    export_triangle_info(points, triangles, neighbors, out_file="mesh_info.txt")
    export_csv(points, triangles, neighbors, out_file="mesh_info.csv")

    print("Fichiers générés :")
    print("  - mesh_info.txt")
    print("  - mesh_info.csv")


if __name__ == "__main__":
    main()