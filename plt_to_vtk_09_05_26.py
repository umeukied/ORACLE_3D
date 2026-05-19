#!/usr/bin/env python3
"""
Convertit field.plt (Tecplot ASCII multi-zones) en serie VTK pour ParaView.
Tous les pas de temps sont convertis et places dans un sous-dossier.

Usage : python3 plt_to_vtk.py field.plt
Resultat :
  - field_vtk/field_0000.vts
  - field_vtk/field_0001.vts
  - ...
  - field_vtk/field.pvd  <=  ouvrir CE fichier dans ParaView
"""

import sys
import re
import os
import numpy as np


def parse_plt(plt_file):
    variables = []
    zones = []

    with open(plt_file, 'r') as f:
        lines = f.readlines()

    # Lire les variables
    for line in lines:
        if 'VARIABLES' in line.upper():
            variables = re.findall(r'"([^"]+)"', line)
            break

    print(f"Variables : {variables}")
    n_vars = len(variables)

    # Lire toutes les zones
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.upper().startswith('ZONE'):
            ni = nj = nk = 1
            time_val = float(len(zones))

            m = re.search(r'T="time=\s*([\d.eE+\-]+)"', line)
            if m: time_val = float(m.group(1))
            m = re.search(r'\bI\s*=\s*(\d+)', line, re.IGNORECASE)
            if m: ni = int(m.group(1))
            m = re.search(r'\bJ\s*=\s*(\d+)', line, re.IGNORECASE)
            if m: nj = int(m.group(1))
            m = re.search(r'\bK\s*=\s*(\d+)', line, re.IGNORECASE)
            if m: nk = int(m.group(1))

            n_pts = ni * nj * nk
            i += 1

            raw = []
            while i < len(lines) and len(raw) < n_pts:
                vals = lines[i].split()
                if vals:
                    try:
                        row = [float(v) for v in vals]
                        if len(row) == n_vars:
                            raw.append(row)
                    except ValueError:
                        pass
                i += 1

            if len(raw) == n_pts:
                zones.append({
                    'ni': ni, 'nj': nj, 'nk': nk,
                    'time': time_val,
                    'data': np.array(raw)
                })
                print(f"  Zone {len(zones):3d} : {ni}x{nj}x{nk}  t={time_val:.4e}")
            else:
                print(f"  Zone incomplete ({len(raw)}/{n_pts}), ignoree")
        else:
            i += 1

    return variables, zones


def write_vts(filepath, zone, variables):
    ni, nj, nk = zone['ni'], zone['nj'], zone['nk']
    data = zone['data']
    n_pts = ni * nj * nk

    var_lower = [v.lower() for v in variables]
    u_idx = next((i for i,v in enumerate(var_lower) if v=='u'), None)
    v_idx = next((i for i,v in enumerate(var_lower) if v=='v'), None)
    w_idx = next((i for i,v in enumerate(var_lower) if v=='w'), None)
    p_idx = next((i for i,v in enumerate(var_lower) if v=='p'), None)

    with open(filepath, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write(f'  <StructuredGrid WholeExtent="0 {ni-1} 0 {nj-1} 0 {nk-1}">\n')
        f.write(f'    <Piece Extent="0 {ni-1} 0 {nj-1} 0 {nk-1}">\n')

        # Coordonnees
        f.write('      <Points>\n')
        f.write('        <DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for pt in range(n_pts):
            f.write(f'          {data[pt,0]:.6e} {data[pt,1]:.6e} {data[pt,2]:.6e}\n')
        f.write('        </DataArray>\n')
        f.write('      </Points>\n')

        f.write('      <PointData>\n')

        # Vecteur vitesse
        if all(x is not None for x in [u_idx, v_idx, w_idx]):
            f.write('        <DataArray type="Float32" Name="Velocity" '
                    'NumberOfComponents="3" format="ascii">\n')
            for pt in range(n_pts):
                f.write(f'          {data[pt,u_idx]:.6e} {data[pt,v_idx]:.6e} {data[pt,w_idx]:.6e}\n')
            f.write('        </DataArray>\n')

        # Pression
        if p_idx is not None:
            f.write('        <DataArray type="Float32" Name="Pressure" format="ascii">\n')
            for pt in range(n_pts):
                f.write(f'          {data[pt,p_idx]:.6e}\n')
            f.write('        </DataArray>\n')

        # Autres scalaires
        skip = {0, 1, 2}
        for idx in [u_idx, v_idx, w_idx, p_idx]:
            if idx is not None:
                skip.add(idx)
        for col_idx, vname in enumerate(variables):
            if col_idx in skip:
                continue
            f.write(f'        <DataArray type="Float32" Name="{vname}" format="ascii">\n')
            for pt in range(n_pts):
                f.write(f'          {data[pt,col_idx]:.6e}\n')
            f.write('        </DataArray>\n')

        f.write('      </PointData>\n')
        f.write('    </Piece>\n')
        f.write('  </StructuredGrid>\n')
        f.write('</VTKFile>\n')


def write_pvd(pvd_path, vts_names, times):
    with open(pvd_path, 'w') as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1">\n')
        f.write('  <Collection>\n')
        for name, t in zip(vts_names, times):
            f.write(f'    <DataSet timestep="{t:.6e}" file="{name}"/>\n')
        f.write('  </Collection>\n')
        f.write('</VTKFile>\n')


def convert(plt_file):
    # Creer le dossier de sortie : field_vtk/ (ou vortex_vtk/ etc.)
    base      = os.path.splitext(os.path.basename(plt_file))[0]
    out_dir   = os.path.join(os.path.dirname(plt_file), base + '_vtk')
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n=== Dossier de sortie : {out_dir} ===\n")

    variables, zones = parse_plt(plt_file)
    if not zones:
        print("Aucune zone trouvee.")
        return

    vts_names = []
    times     = []

    for idx, zone in enumerate(zones):
        vts_name = f"{base}_{idx:04d}.vts"
        vts_path = os.path.join(out_dir, vts_name)
        write_vts(vts_path, zone, variables)
        vts_names.append(vts_name)
        times.append(zone['time'])

    pvd_path = os.path.join(out_dir, base + '.pvd')
    write_pvd(pvd_path, vts_names, times)

    print(f"\n=== Termine ===")
    print(f"  {len(zones)} fichiers .vts dans : {out_dir}/")
    print(f"  Fichier collection  : {pvd_path}")
    print(f"\n=> Lance ParaView avec :")
    print(f"   paraview {pvd_path}")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage : python3 plt_to_vtk.py field.plt")
        sys.exit(1)
    convert(sys.argv[1])