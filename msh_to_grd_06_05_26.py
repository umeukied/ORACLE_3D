#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
msh_to_oracle_grd_mendez2020.py — version corrigee

Correction : ajout de "symmetry" et "symmetry_inner" dans les
candidats pour le groupe 'slip' (paroi interieure type 3).
"""

import argparse
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

import numpy as np


class BlockDims:
    def __init__(self, ni, nj, nk):
        self.ni = int(ni)
        self.nj = int(nj)
        self.nk = int(nk)

    @property
    def nim(self): return self.ni + 1
    @property
    def njm(self): return self.nj + 1
    @property
    def nkm(self): return self.nk + 1
    @property
    def ncv(self): return self.ni * self.nj * self.nk


def oracle_inp_3d(i, j, k, dims):
    return (j + 1) + (i) * dims.nje + (k) * dims.nie * dims.nje

@property
def nje(self): return self.nj + 2
@property
def nie(self): return self.ni + 2

BlockDims.nje = property(lambda self: self.nj + 2)
BlockDims.nie = property(lambda self: self.ni + 2)


HEX8_FACES = [
    (0, 1, 2, 3), (4, 5, 6, 7),
    (0, 1, 5, 4), (1, 2, 6, 5),
    (2, 3, 7, 6), (3, 0, 4, 7),
]


def find_section(lines, name):
    start = end = None
    for i, line in enumerate(lines):
        if line.strip() == f"${name}":     start = i
        elif line.strip() == f"$End{name}": end = i; break
    return start, end


def parse_mesh_format(lines):
    s, _ = find_section(lines, "MeshFormat")
    if s is None: raise ValueError("Missing $MeshFormat")
    parts = lines[s + 1].split()
    if int(parts[1]) != 0: raise ValueError("ASCII only")
    return parts[0]


def parse_physical_names(lines):
    physical_names = {}
    s, _ = find_section(lines, "PhysicalNames")
    if s is None: return physical_names
    n = int(lines[s + 1].strip())
    for k in range(n):
        parts = lines[s + 2 + k].split()
        physical_names[int(parts[1])] = " ".join(parts[2:]).strip().strip('"')
    return physical_names


def read_msh2_nodes(lines):
    s, _ = find_section(lines, "Nodes")
    if s is None: s, _ = find_section(lines, "ParametricNodes")
    if s is None: raise ValueError("Missing $Nodes")
    n = int(lines[s + 1].strip())
    nodes = {}
    for k in range(n):
        p = lines[s + 2 + k].split()
        nodes[int(p[0])] = np.array(list(map(float, p[1:4])))
    return nodes


def read_msh2_ascii(lines):
    nodes = read_msh2_nodes(lines)
    s, _ = find_section(lines, "Elements")
    if s is None: raise ValueError("Missing $Elements")
    n = int(lines[s + 1].strip())
    hexes, quads = [], []
    for k in range(n):
        p = lines[s + 2 + k].split()
        eid, etype, ntags = int(p[0]), int(p[1]), int(p[2])
        tags = list(map(int, p[3:3 + ntags]))
        conn = list(map(int, p[3 + ntags:]))
        phys = tags[0] if ntags >= 1 else 0
        if   etype == 3 and len(conn) == 4: quads.append({"id": eid, "phys": phys, "nodes": conn})
        elif etype == 5 and len(conn) == 8: hexes.append({"id": eid, "phys": phys, "nodes": conn})
    print("Diagnostic MSH2:")
    print("  Physical tags quads :", dict(Counter(q["phys"] for q in quads)))
    print("  Physical tags hexes :", dict(Counter(h["phys"] for h in hexes)))
    return nodes, hexes, quads


def parse_entities_msh4(lines):
    mapping = {}
    s, _ = find_section(lines, "Entities")
    if s is None: return mapping
    header = list(map(int, lines[s + 1].split()))
    num_points, num_curves, num_surfaces, num_volumes = header
    idx = s + 2
    for dim, count in [(0, num_points), (1, num_curves), (2, num_surfaces), (3, num_volumes)]:
        for _ in range(count):
            parts = lines[idx].split(); idx += 1
            tag = int(parts[0])
            np_pos = 4 if dim == 0 else 7
            num_phys = int(parts[np_pos]) if len(parts) > np_pos else 0
            phys = list(map(int, parts[np_pos + 1:np_pos + 1 + num_phys]))
            mapping[(dim, tag)] = phys
    return mapping


def read_msh4_ascii(lines):
    s, _ = find_section(lines, "Nodes")
    if s is None: raise ValueError("Missing $Nodes")
    header = lines[s + 1].split()
    num_entity_blocks = int(header[0])
    total_num_nodes   = int(header[1])
    nodes = {}
    idx = s + 2
    for _ in range(num_entity_blocks):
        parts = lines[idx].split(); num_in_block = int(parts[3]); idx += 1
        tags = [int(lines[idx + i].strip()) for i in range(num_in_block)]; idx += num_in_block
        for local in range(num_in_block):
            vals = list(map(float, lines[idx].split())); idx += 1
            nodes[tags[local]] = np.array(vals[:3])
    s, _ = find_section(lines, "Elements")
    if s is None: raise ValueError("Missing $Elements")
    num_entity_blocks = int(lines[s + 1].split()[0])
    hexes, quads = [], []
    idx = s + 2
    for _ in range(num_entity_blocks):
        parts = lines[idx].split()
        entity_dim, entity_tag, element_type, num_in_block = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
        idx += 1
        for _ in range(num_in_block):
            vals = list(map(int, lines[idx].split())); idx += 1
            eid, conn = vals[0], vals[1:]
            if   element_type == 3 and len(conn) == 4:
                quads.append({"id": eid, "phys": entity_tag, "nodes": conn, "entity_tag": entity_tag})
            elif element_type == 5 and len(conn) == 8:
                hexes.append({"id": eid, "phys": entity_tag, "nodes": conn, "entity_tag": entity_tag})
    return nodes, hexes, quads


def read_msh_ascii(msh_path):
    with open(msh_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [line.strip() for line in f]
    version        = parse_mesh_format(lines)
    physical_names = parse_physical_names(lines)
    if version.startswith("2."):
        nodes, hexes, quads = read_msh2_ascii(lines)
    elif version.startswith("4."):
        entity_phys_map = parse_entities_msh4(lines)
        nodes, hexes, quads = read_msh4_ascii(lines)
        for q in quads:
            pl = entity_phys_map.get((2, q["entity_tag"]), [])
            q["phys"] = pl[0] if pl else q["phys"]
        for h in hexes:
            pl = entity_phys_map.get((3, h["entity_tag"]), [])
            h["phys"] = pl[0] if pl else h["phys"]
    else:
        raise ValueError(f"Version non supportee: {version}")
    if not hexes: raise ValueError("Aucun Hexa8 trouve.")
    if not quads: raise ValueError("Aucun Quad4 trouve.")
    return version, physical_names, nodes, hexes, quads


def build_face_map_and_adjacency(hexes):
    face_map  = defaultdict(list)
    for cidx, h in enumerate(hexes):
        hn = h["nodes"]
        for lf, loc in enumerate(HEX8_FACES):
            key = tuple(sorted(hn[i] for i in loc))
            face_map[key].append((cidx, lf))
    adjacency = defaultdict(list)
    for owners in face_map.values():
        if len(owners) == 2:
            c0, c1 = owners[0][0], owners[1][0]
            adjacency[c0].append(c1); adjacency[c1].append(c0)
    return face_map, adjacency


def quad_centroid(nodes, nids):
    return np.mean([nodes[n] for n in nids], axis=0)


def quad_normal_and_area(nodes, nids):
    p0, p1, p2, p3 = [np.array(nodes[n]) for n in nids]
    n1 = np.cross(p1 - p0, p2 - p0); n2 = np.cross(p2 - p0, p3 - p0)
    n = n1 + n2; area = 0.5 * np.linalg.norm(n1) + 0.5 * np.linalg.norm(n2)
    nn = np.linalg.norm(n)
    return (n / nn if nn > 1e-30 else np.zeros(3)), area


def build_boundary_items(quads, nodes, face_map, physical_names):
    items = []
    for q in quads:
        fkey = tuple(sorted(q["nodes"]))
        owners = face_map.get(fkey, [])
        if len(owners) != 1: continue
        cidx = owners[0][0]
        phys_name = physical_names.get(q["phys"], f"phys_{q['phys']}")
        items.append({"quad": q, "cell": cidx,
                      "centroid": quad_centroid(nodes, q["nodes"]),
                      "normal": quad_normal_and_area(nodes, q["nodes"])[0],
                      "area": quad_normal_and_area(nodes, q["nodes"])[1],
                      "phys_name": phys_name})
    if not items: raise ValueError("Aucun quad de bord valide.")
    return items


def normalize_name(name):
    return str(name).strip().lower().replace("-", "_").replace(" ", "_")


def first_non_empty(groups, candidates):
    for c in candidates:
        if groups.get(c): return groups[c]
    return []


def classify_boundary_quads(quads, physical_names, nodes, face_map):
    items = build_boundary_items(quads, nodes, face_map, physical_names)
    groups_raw = defaultdict(list)
    for it in items:
        groups_raw[normalize_name(it["phys_name"])].append(it)

    print("Groupes physiques detectes :")
    for name, vals in sorted(groups_raw.items()):
        print(f"  {name:35s} : {len(vals)}")

    inlet_cands  = ["inlet", "in", "entree"]
    outlet_cands = ["outlet", "out", "sortie"]

    # j=1 : paroi concave no-slip
    wall_cands = [
        "wall_concave", "concave", "wall_bottom", "bottom",
        "lower_wall", "wall_lower", "wall"
    ]

    # -------------------------------------------------------
    # CORRECTION : ajout de "symmetry", "symmetry_inner",
    # "symmetry_top" dans les candidats pour la paroi
    # interieure (anciennement wall_slip)
    # -------------------------------------------------------
    slip_cands = [
        "wall_slip", "slip",
        "symmetry",          # <-- NOUVEAU : nom utilise dans le .geo
        "symmetry_inner",    # <-- NOUVEAU : variante explicite
        "symmetry_top",      # <-- NOUVEAU : variante
        "inner",
        "wall_top", "top_wall", "upper_wall", "top",
        "farfield_top",
    ]

    zmin_cands = [
        "symmetry_z0", "symmetry_zmin", "zmin", "side1",
        "wall_lateral_1", "lateral_1", "symmetry1"
    ]
    zmax_cands = [
        "symmetry_zl", "symmetry_zlz", "symmetry_zmax", "zmax", "side2",
        "wall_lateral_2", "lateral_2", "symmetry2"
    ]

    groups = defaultdict(list)
    groups["inlet"]         = first_non_empty(groups_raw, inlet_cands)
    groups["outlet"]        = first_non_empty(groups_raw, outlet_cands)
    groups["wall"]          = first_non_empty(groups_raw, wall_cands)
    groups["slip"]          = first_non_empty(groups_raw, slip_cands)
    groups["symmetry_zmin"] = first_non_empty(groups_raw, zmin_cands)
    groups["symmetry_zmax"] = first_non_empty(groups_raw, zmax_cands)

    required = ["inlet", "outlet", "wall", "slip", "symmetry_zmin", "symmetry_zmax"]
    for name in required:
        if not groups[name]:
            raise ValueError(
                f"Groupe '{name}' vide.\n"
                f"Groupes detectes : {list(groups_raw.keys())}\n"
                f"Candidats testes : voir script."
            )
    return groups


def check_boundary_group_consistency(boundary_groups):
    nin   = len(boundary_groups["inlet"])
    nout  = len(boundary_groups["outlet"])
    nwall = len(boundary_groups["wall"])
    nslip = len(boundary_groups["slip"])
    nzmin = len(boundary_groups["symmetry_zmin"])
    nzmax = len(boundary_groups["symmetry_zmax"])
    print(f"  inlet/outlet          : {nin}/{nout}")
    print(f"  wall_concave/slip     : {nwall}/{nslip}")
    print(f"  symmetry_z0/zL        : {nzmin}/{nzmax}")
    if nin != nout:  raise ValueError(f"inlet ({nin}) != outlet ({nout})")
    if nwall != nslip: raise ValueError(f"wall ({nwall}) != slip ({nslip})")
    if nzmin != nzmax: raise ValueError(f"z0 ({nzmin}) != zL ({nzmax})")


def bfs_distances(adjacency, sources, ncells):
    dist = np.full(ncells, -1, dtype=int)
    q = deque()
    for s in sources:
        dist[s] = 0; q.append(s)
    while q:
        c = q.popleft()
        for nb in adjacency[c]:
            if dist[nb] == -1:
                dist[nb] = dist[c] + 1; q.append(nb)
    return dist


def reconstruct_cell_indices(hexes, adjacency, boundary_groups):
    ncells = len(hexes)
    def cells(g): return sorted(set(it["cell"] for it in boundary_groups[g]))
    di = bfs_distances(adjacency, cells("inlet"),         ncells)
    dj = bfs_distances(adjacency, cells("wall"),          ncells)
    dk = bfs_distances(adjacency, cells("symmetry_zmin"), ncells)
    if np.any(di < 0) or np.any(dj < 0) or np.any(dk < 0):
        raise ValueError("Cellules inaccessibles depuis un bord.")
    i_idx = di + 1; j_idx = dj + 1; k_idx = dk + 1
    ni = int(np.max(i_idx)); nj = int(np.max(j_idx)); nk = int(np.max(k_idx))
    dims = BlockDims(ni, nj, nk)
    # Verifications
    for name, expected, arr in [
        ("inlet",         1,       i_idx),
        ("outlet",        ni,      i_idx),
        ("wall",          1,       j_idx),
        ("slip",          nj,      j_idx),
        ("symmetry_zmin", 1,       k_idx),
        ("symmetry_zmax", nk,      k_idx),
    ]:
        if any(arr[c] != expected for c in cells(name)):
            raise ValueError(f"Coherence echouee pour '{name}'.")
    cell_logical = {}; seen = set()
    for c in range(ncells):
        key = (int(i_idx[c]), int(j_idx[c]), int(k_idx[c]))
        if key in seen: raise ValueError(f"Doublon : {key}")
        seen.add(key); cell_logical[c] = key
    if len(seen) != dims.ncv:
        raise ValueError(f"Reconstruction incomplete : {len(seen)} != {dims.ncv}")
    return dims, cell_logical


def reconstruct_vertex_grid(nodes, hexes, boundary_groups, dims, cell_logical):
    def nset(g): return set(n for it in boundary_groups[g] for n in it["quad"]["nodes"])
    inlet_n = nset("inlet"); outlet_n = nset("outlet")
    wall_n  = nset("wall");  slip_n   = nset("slip")
    zmin_n  = nset("symmetry_zmin"); zmax_n = nset("symmetry_zmax")
    node_adj = defaultdict(list)
    for cidx, h in enumerate(hexes):
        li, lj, lk = cell_logical[cidx]
        for nid in h["nodes"]: node_adj[nid].append((li, lj, lk))
    V = np.full((dims.nkm, dims.nim, dims.njm, 3), np.nan)
    for nid, coord in nodes.items():
        if nid not in node_adj: continue
        adj = node_adj[nid]
        is_ = [a[0] for a in adj]; js_ = [a[1] for a in adj]; ks_ = [a[2] for a in adj]
        iv = 1 if nid in inlet_n  else (dims.ni+1 if nid in outlet_n else max(is_))
        jv = 1 if nid in wall_n   else (dims.nj+1 if nid in slip_n   else max(js_))
        kv = 1 if nid in zmin_n   else (dims.nk+1 if nid in zmax_n   else max(ks_))
        kk, ii, jj = kv-1, iv-1, jv-1
        if not np.isnan(V[kk, ii, jj, 0]):
            if np.linalg.norm(V[kk, ii, jj] - coord) > 1e-12:
                raise ValueError(f"Conflit sommet ({kv},{iv},{jv})")
        else:
            V[kk, ii, jj] = coord
    missing = np.argwhere(np.isnan(V[..., 0]))
    if len(missing): raise ValueError(f"{len(missing)} sommets manquants.")
    return V


def build_boundary_regions(dims):
    regions = []
    # IDR 1 : inlet (type 1, ISIDE=1)
    faces = [(oracle_inp_3d(1, j, k, dims), 1)
             for k in range(1, dims.nk+1) for j in range(1, dims.nj+1)]
    regions.append((1, 1, faces))
    # IDR 2 : outlet (type 2, ISIDE=2)
    faces = [(oracle_inp_3d(dims.ni, j, k, dims), 2)
             for k in range(1, dims.nk+1) for j in range(1, dims.nj+1)]
    regions.append((2, 2, faces))
    # IDR 3 : wall_concave (type 4, ISIDE=3, j=1)
    faces = [(oracle_inp_3d(i, 1, k, dims), 3)
             for k in range(1, dims.nk+1) for i in range(1, dims.ni+1)]
    regions.append((3, 4, faces))
    # IDR 4 : symmetry paroi interne (type 3, ISIDE=4, j=nj)
    faces = [(oracle_inp_3d(i, dims.nj, k, dims), 4)
             for k in range(1, dims.nk+1) for i in range(1, dims.ni+1)]
    regions.append((4, 3, faces))
    # IDR 5 : symmetry_z0 (type 3, ISIDE=5, k=1)
    faces = [(oracle_inp_3d(i, j, 1, dims), 5)
             for j in range(1, dims.nj+1) for i in range(1, dims.ni+1)]
    regions.append((5, 3, faces))
    # IDR 6 : symmetry_zL (type 3, ISIDE=6, k=nk)
    faces = [(oracle_inp_3d(i, j, dims.nk, dims), 6)
             for j in range(1, dims.nj+1) for i in range(1, dims.ni+1)]
    regions.append((6, 3, faces))
    return regions


def write_oracle_grd(out_path, V, dims, boundary_regions):
    with open(out_path, "w") as f:
        f.write("1\n")
        f.write(f"{dims.ncv}\n")
        f.write("1\n")
        f.write(f"{dims.ni} {dims.nj} {dims.nk}\n")
        f.write("0\n")
        f.write(f"{len(boundary_regions)}\n")
        for kk in range(dims.nkm):
            for ii in range(dims.nim):
                for jj in range(dims.njm):
                    x, y, z = V[kk, ii, jj]
                    f.write(f"{x:.16e} {y:.16e} {z:.16e}\n")
        for rid, bc_type, faces in boundary_regions:
            f.write(f"{rid} {bc_type} {len(faces)}\n")
            for inp, iside in faces:
                f.write(f"{inp} {iside}\n")
    print(f"  {dims.ncv} volumes de controle")
    print(f"  {dims.nim}x{dims.njm}x{dims.nkm} sommets")
    print(f"  {len(boundary_regions)} regions de bord")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("msh", type=Path)
    parser.add_argument("grd", type=Path)
    args = parser.parse_args()
    if not args.msh.exists(): raise FileNotFoundError(f"Introuvable : {args.msh}")

    print("="*60)
    print("Lecture du .msh ...")
    version, physical_names, nodes, hexes, quads = read_msh_ascii(args.msh)
    print(f"Version : {version}")
    print(f"Noms physiques : {physical_names}")
    print(f"Noeuds : {len(nodes)}  Hexaedres : {len(hexes)}  Quads : {len(quads)}")

    face_map, adjacency = build_face_map_and_adjacency(hexes)

    print("Classification des frontieres ...")
    boundary_groups = classify_boundary_quads(quads, physical_names, nodes, face_map)
    check_boundary_group_consistency(boundary_groups)

    print("Reconstruction logique ...")
    dims, cell_logical = reconstruct_cell_indices(hexes, adjacency, boundary_groups)
    print(f"Bloc : ni={dims.ni}  nj={dims.nj}  nk={dims.nk}")

    V = reconstruct_vertex_grid(nodes, hexes, boundary_groups, dims, cell_logical)
    boundary_regions = build_boundary_regions(dims)

    print("Ecriture .grd ...")
    write_oracle_grd(args.grd, V, dims, boundary_regions)

    print("="*60)
    print("Succes !")
    print(f"  Entree  : {args.msh}")
    print(f"  Sortie  : {args.grd}")
    print("Regions Oracle3D :")
    print("  IDR 1  inlet        -> type 1  INLET")
    print("  IDR 2  outlet       -> type 2  OUTLET")
    print("  IDR 3  wall_concave -> type 4  WALL")
    print("  IDR 4  symmetry     -> type 3  SYMMETRY (paroi interne)")
    print("  IDR 5  symmetry_z0  -> type 3  SYMMETRY")
    print("  IDR 6  symmetry_zL  -> type 3  SYMMETRY")
    print("="*60)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[ERREUR] {e}", file=sys.stderr)
        sys.exit(1)
