///////////////////////////////////////////////////////////////
// gortler_top_inlet_bottom_outlet.geo
//
// Canal courbe de Gortler - Conditions aux limites modifiees :
//
//   inlet        -> type 1  (entree, face horizontale en haut)
//   outlet       -> type 2  (sortie, face verticale a droite)
//   wall_concave -> type 4  (paroi solide, arc exterieur r=R+H)
//   symmetry     -> type 3  (paroi interieure, arc interieur r=R)
//   symmetry_z0  -> type 3  (face laterale z=0)
//   symmetry_zL  -> type 3  (face laterale z=Lz)
//
// Geometrie :
//   Centre du coude C = (R+H, R+H)
//   Arc exterieur r=R+H : paroi concave (wall_concave, type 4)
//   Arc interieur r=R   : symetrie      (symmetry,    type 3)
//   Inlet  = segment horizontal en haut (y=R+H, x: 0..H)
//   Outlet = segment vertical a droite  (x=R+H, y: 0..H)
///////////////////////////////////////////////////////////////

SetFactory("OpenCASCADE");

// ------------------------------------------------------------
// Parametres geometriques
// ------------------------------------------------------------
R   = 30.0;   // rayon interieur (paroi symetrie)
H   = 10.0;   // epaisseur du canal
Lz  = 20.0;   // profondeur spanwise

// ------------------------------------------------------------
// Parametres de maillage
// ------------------------------------------------------------
n_eta = 50;       // direction normale wall -> symmetry
bump_eta = 0.0001; // concentration forte vers paroi concave (j=1)
n_xi  = 120;       // direction streamwise
n_z   = 70;       // direction spanwise

lc = 1.0;

// ------------------------------------------------------------
// Geometrie
// Centre en (R+H, R+H) => inlet en haut, outlet a droite
// ------------------------------------------------------------
Cx = R + H;
Cy = R + H;

p_center = newp; Point(p_center) = {Cx, Cy, 0, lc};

// Inlet (haut, y = R+H)
p_in_wall = newp; Point(p_in_wall) = {0,   R+H, 0, lc}; // wall_concave @ inlet
p_in_slip = newp; Point(p_in_slip) = {H,   R+H, 0, lc}; // symmetry    @ inlet

// Outlet (droite, x = R+H)
p_out_slip = newp; Point(p_out_slip) = {R+H, H, 0, lc}; // symmetry    @ outlet
p_out_wall = newp; Point(p_out_wall) = {R+H, 0, 0, lc}; // wall_concave @ outlet

// ------------------------------------------------------------
// Courbes
// ------------------------------------------------------------

// Inlet horizontal en haut
c_inlet = newl; Line(c_inlet) = {p_in_wall, p_in_slip};

// Arc interieur r=R : symetrie
c_arc_slip = newl; Circle(c_arc_slip) = {p_in_slip, p_center, p_out_slip};

// Outlet vertical a droite
c_outlet = newl; Line(c_outlet) = {p_out_slip, p_out_wall};

// Arc exterieur r=R+H : wall_concave
c_arc_wall = newl; Circle(c_arc_wall) = {p_out_wall, p_center, p_in_wall};

// ------------------------------------------------------------
// Surface 2D
// ------------------------------------------------------------
Curve Loop(1) = {c_inlet, c_arc_slip, c_outlet, c_arc_wall};
Plane Surface(1) = {1};

// ------------------------------------------------------------
// Maillage structure transfinite
// ------------------------------------------------------------

// Streamwise : uniforme sur les deux arcs
Transfinite Curve {c_arc_wall, c_arc_slip} = n_xi + 1 Using Bump 1.0;

// Wall-normal : concentration forte vers p_in_wall (paroi concave, j=1)
// c_inlet : de p_in_wall (wall) vers p_in_slip (symmetry)
// bump_eta tres petit => premiere cellule tres fine cote wall
Transfinite Curve {c_inlet, c_outlet} = n_eta + 1 Using Bump bump_eta;

// Coins : wall_concave inlet -> symmetry inlet -> symmetry outlet -> wall_concave outlet
Transfinite Surface {1} = {p_in_wall, p_in_slip, p_out_slip, p_out_wall};
Recombine Surface {1};

// ------------------------------------------------------------
// Extrusion 3D selon z
// ------------------------------------------------------------
out[] = Extrude {0, 0, Lz} {
  Surface{1};
  Layers{n_z};
  Recombine;
};

// out[0] = face z=Lz (symmetry_zL)
// out[1] = volume fluide
// out[2] = extrusion de c_inlet      => inlet
// out[3] = extrusion de c_arc_slip   => symmetry (arc interieur r=R)
// out[4] = extrusion de c_outlet     => outlet
// out[5] = extrusion de c_arc_wall   => wall_concave (arc exterieur r=R+H)

// ------------------------------------------------------------
// Options de maillage
// ------------------------------------------------------------
Mesh.Algorithm    = 8;
Mesh.Algorithm3D  = 1;
Mesh.RecombineAll           = 1;
Mesh.RecombinationAlgorithm = 1;
Mesh.Optimize       = 0;
Mesh.OptimizeNetgen = 0;
Mesh.Smoothing      = 0;
Mesh.MshFileVersion = 2.2;
Mesh.SaveAll        = 0;

// ------------------------------------------------------------
// Groupes physiques
//
// Nommage reconnu par le convertisseur Oracle3D :
//   inlet        -> type 1
//   outlet       -> type 2
//   symmetry_z0  -> type 3
//   symmetry_zL  -> type 3
//   symmetry     -> type 3  (arc interieur r=R)
//   wall_concave -> type 4  (arc exterieur r=R+H)
// ------------------------------------------------------------
Physical Volume("fluid") = {out[1]};

Physical Surface("inlet")        = {out[2]};  // type 1
Physical Surface("symmetry")     = {out[3]};  // type 3 - arc interieur r=R
Physical Surface("outlet")       = {out[4]};  // type 2
Physical Surface("wall_concave") = {out[5]};  // type 4 - arc exterieur r=R+H
Physical Surface("symmetry_z0")  = {1};       // type 3
Physical Surface("symmetry_zL")  = {out[0]}; // type 3

///////////////////////////////////////////////////////////////
// RESUME
//
// Paroi exterieure (r=R+H=40) : wall_concave -> type 4 (no-slip)
// Paroi interieure (r=R=30)   : symmetry     -> type 3
// Inlet  (y=R+H, horizontal)  : inlet        -> type 1
// Outlet (x=R+H, vertical)    : outlet       -> type 2
// Faces spanwise              : symmetry_z0/zL -> type 3
//
// Maillage :
//   n_eta = 120, n_xi = 300, n_z = 200
//   Total = 7 200 000 cellules
//   Premiere cellule : ~0.005 delta0 (bump_eta=0.0001)
//
// Commande :
//   gmsh gortler_top_inlet_bottom_outlet.geo -3 -format msh2 \
//        -o gortler.msh
///////////////////////////////////////////////////////////////