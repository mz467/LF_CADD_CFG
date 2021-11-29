# Coupled LAMMPS and FEniCS
# To convert the atoms information generated from mesh.py to LAMMPS data format
# Written by WenjiaG. on 08/13/2018

import math
import os
from clf_kfield_displacement import kfield_disp
import sys
from clf_key_word_search import *

# ====================== Initial atom data preparation for LAMMPS ============================== #

def lammps_intial_atom_data(atom_filename, crack_x, crack_y, K_I, E, nu, kappa):
    interface_coord = []
    pad_atom_coord = []

    # Check the number and the position of atoms
    xmax = 0
    xmin = 0
    ymax = 0
    ymin = 0

    f = open(atom_filename, 'r')

    fw = open('atoms.data', 'w')

    for n, line in enumerate(f):
        line = line.strip()
        columns = line.split()
        atom_id = int(columns[0])
        atom_type = int(columns[1])
        atom_x = float(columns[2])
        atom_y = float(columns[3])
        atom_z = float(columns[4])
        atom_a = columns[5]
        atom_b = columns[6]
        atom_c = columns[7]

        if atom_type == 3:
            pad_atom_coord.append([int(atom_id), float(atom_x), float(atom_y)])

        elif atom_type == 2:
            interface_coord.append([int(atom_id), float(atom_x), float(atom_y)])

        # --- apply linear elastic displacement field to all atoms ----
        u_x, u_y = kfield_disp(atom_x, atom_y, crack_x, crack_y, K_I, E, nu, kappa)
        atom_x = atom_x + u_x
        atom_y = atom_y + u_y

        fw.write('%s %s %s %s %s %s %s %s\n' % (str(atom_id), str(atom_type), str(atom_x), str(atom_y),
                                                str(atom_z), atom_a, atom_b, atom_c))

        if atom_x > xmax:
            xmax = atom_x
        elif atom_x < xmin:
            xmin = atom_x
        if atom_y > ymax:
            ymax = atom_y
        elif atom_y < ymin:
            ymin = atom_y

    n_atoms = n + 1
    f.close()
    fw.close()

    # Calculate xlo, xhi, ylo, yhi
    xlo = xmin - (xmax - xmin) / 10000
    xhi = xmax + (xmax - xmin) / 10000
    ylo = ymin - (ymax - ymin) / 10000
    yhi = ymax + (ymax - ymin) / 10000

    # Copy the content to a new file and add more information required by LAMMPS
    f2 = open("atoms.data", "r")
    contents = f2.readlines()
    f2.close()

    contents.insert(0, 'LAMMPS data file converted from meshing output, timestep = 0 \n')
    contents.insert(1, '\n')
    contents.insert(2, '%s atoms \n' % str(n_atoms))
    contents.insert(3, '3 atom types\n')
    contents.insert(4, '\n')
    contents.insert(5, '%s %s xlo xhi\n' % (str(xlo), str(xhi)))
    contents.insert(6, '%s %s ylo yhi\n' % (str(ylo), str(yhi)))
    contents.insert(7, '-2.7856927797769487e-01 2.7856927797769487e-01 zlo zhi \n')
    contents.insert(8, '\n')
    contents.insert(9, 'Masses \n')
    contents.insert(10, '\n')
    contents.insert(11, '1 1 \n')
    contents.insert(12, '2 1 \n')
    contents.insert(13, '3 1 \n')
    contents.insert(14, '\n')
    contents.insert(15, 'Atoms # atomic \n')
    contents.insert(16, '\n')

    fw2 = open("atoms.lammps", "w")
    contents = "".join(contents)
    fw2.write(contents)
    fw2.close()

    os.remove("atoms.data")

    return interface_coord, n_atoms, pad_atom_coord

# ======================= Update pad atom positions from FEniCS ======================= #

def lammps_pad_update(pad_atom_coord, atoms_data, u_vertex_values):

    n_nodes = int(len(u_vertex_values)/2)

    for i in range(len(pad_atom_coord)):
        atom_id = pad_atom_coord[i][0]
        atom_x = pad_atom_coord[i][1]
        atom_y = pad_atom_coord[i][2]

        u_x = u_vertex_values[atom_id-1]
        u_y = u_vertex_values[atom_id-1+n_nodes]

        atoms_data[3*atom_id-3] = atom_x + u_x
        atoms_data[3*atom_id-2] = atom_y + u_y

    return atoms_data
