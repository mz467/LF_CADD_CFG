# Coupled LAMMPS and FEniCS
# To generate the geometry for mesh generation using Gmsh
# To generate the atoms input data for LAMMPS
# Written by Wenjia G. on 08/29/2018

import math

# ===================== Input parameter ===================== #

# Regions of different sizes of mesh to be generated:
# The first region is the atomistic region, the second one is for pad atoms...
region_x = [150, 160, 180, 220, 300, 460, 780, 1420, 2500]
region_y = [75, 85, 105, 145, 225, 385, 705, 1345, 2500]

# Crack tip position
crack_x_position = -150+25
crack_y_position = 0.0

# Lattice constant
d_equil = 1.1216204059910884 # equilibrium interatomic distance calculated from lammps minimization
d_equil_y = d_equil*math.sqrt(3)/2 # distance between rows of atoms that are parallel to x-axis

crack_x = int(crack_x_position/d_equil)*d_equil
crack_y = int(crack_y_position/d_equil_y)*d_equil_y

# Crack generation:
# Crack tip is at the origin of the sim box
# Three layers of atoms are removed

# ===================== Generate geometry for Gmsh ===================== #
fw = open('300x150-5000_crack-25_mesh.geo', 'w') # output file from this script - the input file for Gmsh
fw1 = open('300x150-5000_crack-25_detection.geo', 'w') # output file from this script - the input file for detection band construction
fw3 = open('300x150-5000_crack-25_atoms.inp', 'w') # atom data file for LAMMPS

id_1234_inner = []
id_5_crack_upper = []
id_9876_outer = []
id_10_outer = []
id_11_crack_lower = []
id_12_inner = []

id_atom_crack_upper = [] # atomistic region meshing - for detection band
id_atom_crack_tip_upper = [] # atomistic region meshing - for detection band
id_atom_crack_tip = [] # atomistic region meshing - for detection band
id_atom_crack_tip_lower = [] # atomistic region meshing - for detection band
id_atom_crack_lower = [] # atomistic region meshing - for detection band

id_atom_point_inside = []

# Boundaries all regions
region_xmin = []
region_xmax = []
region_ymin = []
region_ymax = []

mesh_inside_point_start = []
mesh_inside_point_end = []

node_id = 0 # number of the current point
for region in range(1, len(region_x)+1):

    # =-=-=-=-=-=-=-=-=-=-=-=-= Points generation =-=-=-=-=-=-=-=-=-=-=-=-=
    if region == 1:
        dx = d_equil # distance between two points in x direction in the current region
        dy = d_equil_y # distance between two points in y direction in the current region
    else:
        dx = d_equil*2**(region-2) # distance between two points in x direction in the current region
        dy = d_equil_y*2**(region-2) # distance between two points in x direction in the current region
    nx = int(region_x[region-1]//dx) # number of points to be generated according to the prescribed region size
    ny = int(region_y[region-1]//dy) # number of points to be generated according to the prescribed region size
    if (region == 1 or region == 2) and ny%2 != 0:
        ny = ny+1 # generate even numbers of rows of points for convenience

    # ---------------- Points on boundary --------------
    if region==1 or region==len(region_x):

        # Lower boundary
        node_y = -ny*dy
        for i in range(-nx, nx+1):
            node_x = i*dx
            node_id = node_id + 1
            fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
            if region == 1:
                fw1.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
                fw3.write("%s 2 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y))) # interface nodes
                id_1234_inner.append(node_id)
            elif region == len(region_x):
                id_9876_outer.append(node_id) # outer boundary

        # Right boundary
        for j in range(-ny+1, ny):
            node_y = j * dy
            if region == 1 or region == 2:
                if j%2 == 0:
                    node_x = nx*dx
                else:
                    node_x = dx/2+(nx-1)*dx
            else:
                node_x = nx * dx
            node_id = node_id + 1
            fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
            if region == 1:
                fw1.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
                fw3.write("%s 2 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y))) # interface nodes
                id_1234_inner.append(node_id)
            elif region == len(region_x):
                id_9876_outer.append(node_id) # outer boundary

        # Upper boundary
        node_y = ny*dy
        for i in range(nx, -nx-1, -1):
            node_x = i*dx
            node_id = node_id + 1
            fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
            if region == 1:
                fw1.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
                fw3.write("%s 2 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y))) # interface nodes
                id_1234_inner.append(node_id)
            elif region == len(region_x):
                id_9876_outer.append(node_id) # outer boundary

        ## Crack generation:
        # Left boundary above crack
        for j in range(ny-1, 0, -1):
            node_y = j * dy
            if node_y > 2*d_equil_y:
                if region == 1:
                    if j % 2 == 0:
                        node_x = -nx*dx
                    else:
                        node_x = dx/2-nx*dx
                    node_id = node_id + 1
                    fw1.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
                    fw3.write("%s 2 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y)))
                    id_1234_inner.append(node_id)
                else:
                    node_x = -nx*dx
                    node_id = node_id + 1
                    id_9876_outer.append(node_id)
                fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))

        # Crack upper surface:
        node_y = 2*d_equil_y
        if region == 1:
            node_x = -nx*dx
            node_id = node_id + 1
            fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
            fw1.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
            fw3.write("%s 2 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y)))
            id_5_crack_upper.append(node_id)
        else:
            for i in range(-1,-nx-1,-1):
                node_x = i*dx
                if node_x < region_xmin[region-2]:
                    node_id = node_id+1
                    id_5_crack_upper.append(node_id)
                    fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))

        # Crack lower surface:
        node_y = -2*d_equil_y
        if region == 1:
            node_x = -nx*dx
            node_id = node_id + 1
            fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
            fw1.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
            fw3.write("%s 2 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y)))
            id_11_crack_lower.append(node_id)
        else:
            for i in range(-1,-nx-1,-1):
                node_x = i*dx
                if node_x < region_xmin[region-2]:
                    node_id = node_id+1
                    id_11_crack_lower.append(node_id)
                    fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))

        # Left boundary below crack
        for j in range(0, -ny, -1):
            node_y = j * dy
            if node_y < -2*d_equil_y:
                if region == 1:
                    if j % 2 == 0:
                        node_x = -nx*dx
                    else:
                        node_x = dx/2-nx*dx
                    node_id = node_id + 1
                    fw1.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
                    fw3.write("%s 2 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y)))
                    id_12_inner.append(node_id)
                else:
                    node_x = -nx*dx
                    node_id = node_id + 1
                    id_10_outer.append(node_id)
                fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))


    # Record the boundaries of the current region
    region_xmin.append(-nx * dx)
    region_ymin.append(-ny * dy)
    region_xmax.append(nx * dx)
    region_ymax.append(ny * dy)

    # -------------- Other points generation ----------------
    mesh_inside_point_start.append(node_id+1)

    if region==2:
        for j in range(-ny, ny+1):
            node_y = j * dy
            if j % 2 == 0:
                for i in range(-nx, nx+1):
                    node_x = i * dx
                    # Exclude points that are generated already in the previous regions
                    if not (region_xmin[region-2]<=node_x<= region_xmax[region-2] and region_ymin[region-2]<=node_y<=region_ymax[region-2]):
                        # Exclude points in crack region
                        if not (node_x < 0 and (-2*d_equil_y <= node_y <= 2*d_equil_y)):
                            node_id = node_id + 1
                            fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
                            fw3.write("%s 3 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y)))
            else:
                for i in range(-nx, nx):
                    node_x = dx / 2 + i * dx
                    if not (region_xmin[region-2]<=node_x<=region_xmax[region-2] and region_ymin[region-2]<=node_y<=region_ymax[region-2]):
                        if not (node_x < 0 and (-2 * d_equil_y <= node_y <= 2 * d_equil_y)):
                            node_id = node_id + 1
                            fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
                            fw3.write("%s 3 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y)))
    elif region==len(region_x):
        for j in range(-ny+1, ny):
            node_y = j * dy
            for i in range (-nx+1, nx):
                node_x = i*dx
                if not (region_xmin[region-2]<=node_x<=region_xmax[region-2] and region_ymin[region-2]<=node_y<=region_ymax[region-2]):
                    if not (node_x < 0 and (-2 * d_equil_y <= node_y <= 2 * d_equil_y)):
                        node_id = node_id + 1
                        fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
    elif region != 1:
        for j in range(-ny, ny+1):
            node_y = j * dy
            for i in range (-nx, nx+1):
                node_x = i*dx
                if not (region_xmin[region-2]<=node_x<=region_xmax[region-2] and region_ymin[region-2]<=node_y<=region_ymax[region-2]):
                    if not (node_x < 0 and (-2 * d_equil_y <= node_y <= 2 * d_equil_y)):
                        node_id = node_id + 1
                        fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))

    mesh_inside_point_end.append(node_id)

    if (region != 1 and region != len(region_x)):
        # crack upper surface:
        node_y = 2*d_equil_y
        for i in range(-1, -nx - 1, -1):
            node_x = i*dx
            if node_x < region_xmin[region-2]:
                node_id = node_id+1
                fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
                id_5_crack_upper.append(node_id)
                if region == 2:
                    fw3.write("%s 3 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y)))

        # crack lower surface:
        node_y = -2 * d_equil_y
        for i in range(-1, -nx - 1, -1):
            node_x = i * dx
            if node_x < region_xmin[region-2]:
                node_id = node_id + 1
                fw.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(node_id), str(node_x), str(node_y), str(dx)))
                id_11_crack_lower.append(node_id)
                if region == 2:
                    fw3.write("%s 3 %s %s 0.000000e+00 0 0 0\n" % (str(node_id), str(node_x), str(node_y)))

node_interface = id_1234_inner + [id_5_crack_upper[0]] + [id_11_crack_lower[0]] + id_12_inner
node_outerb = id_9876_outer + [id_5_crack_upper[-1]] + [id_11_crack_lower[-1]] + id_10_outer

# print(id_1234_inner)
# print(id_5_crack_upper)
# print(id_9876_outer)
# print(id_10_outer)
# print(id_11_crack_lower)
# print(id_12_inner)
#print(node_interface)
#print(node_outerb)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Line generation =-=-=-=-=-=-=-=-=-=-=-=-=
phy_line_inner1_min = 1
for i in range (0, len(id_1234_inner)-1):
    line_id = id_1234_inner[i]
    fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_1234_inner[i]), str(id_1234_inner[i+1])))

for i in range (-1, len(id_5_crack_upper)-1):
    line_id = line_id+1
    if i == -1:
        fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_1234_inner[-1]), str(id_5_crack_upper[0])))
        phy_line_inner1_max = line_id
    else:
        fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_5_crack_upper[i]), str(id_5_crack_upper[i+1])))

for i in range(0, len(id_9876_outer)):
    line_id = line_id+1
    if i == 0:
        fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_5_crack_upper[-1]), str(id_9876_outer[-1])))
        phy_line_outer_min = line_id
    else:
        fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_9876_outer[-i]), str(id_9876_outer[-i-1])))

for i in range(0, len(id_10_outer)):
    line_id = line_id+1
    if i == 0:
        fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_9876_outer[0]), str(id_10_outer[-1])))
    else:
        fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_10_outer[-i]), str(id_10_outer[-i-1])))

for i in range(0, len(id_11_crack_lower)):
    line_id = line_id+1
    if i == 0:
        fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_10_outer[0]), str(id_11_crack_lower[-1])))
        phy_line_outer_max = line_id
    else:
        fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_11_crack_lower[-i]), str(id_11_crack_lower[-i-1])))

for i in range(-1, len(id_12_inner)-1):
    line_id = line_id+1
    if i == -1:
        fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_11_crack_lower[0]), str(id_12_inner[0])))
        phy_line_inner2_min = line_id
    else:
        fw.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_12_inner[i]), str(id_12_inner[i+1])))

last_line_id = line_id+1
fw.write("Line(%s) = {%s, %s};\n" % (str(last_line_id), str(id_12_inner[-1]), str(id_1234_inner[0])))
phy_line_inner2_max = last_line_id

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-= Surface generation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

fw.write("Curve Loop(1) = {1:%s};\n" % str(last_line_id))
fw.write("Plane Surface(1) = {1};\n")

# =-=-=-=-=-=-=-=-=-=-=-= Set inside points as part of mesh =-=-=-=-=-=-=-=-=-=-=

for i in range(1, len(mesh_inside_point_start)):
    fw.write("Point {%s:%s} In Surface {1};\n" % (str(mesh_inside_point_start[i]), str(mesh_inside_point_end[i])))

# =-=-=-=-=-=-=-=-=-=-=-=-= Define physical entities =-=-=-=-=-=-=-=-=-=-=-=-=

fw.write("Physical Line(11) = {%s:%s,%s:%s}; \n" % (str(phy_line_inner1_min), str(phy_line_inner1_max),str(phy_line_inner2_min), str(phy_line_inner2_max)))
fw.write("Physical Line(22) = {%s:%s}; \n" % (str(phy_line_outer_min), str(phy_line_outer_max)))

fw.write("Physical Surface(1) = {1}; \n")

fw.close()

# ===================== Generate atoms input data for LAMMPS ===================== #
atom_id = mesh_inside_point_start[2]-1

dx = d_equil  # distance between two points in x direction inside the atomistic region
dy = d_equil_y  # distance between two points in y direction inside the atomistic region

nx = int(region_x[0] // dx)  # number of points to be generated according to the prescribed region size
ny = int(region_y[0] // dy)  # number of points to be generated according to the prescribed region size
if ny % 2 != 0:
    ny = ny + 1  # generate even numbers of rows of points for convenience:

for j in range(-ny + 1, ny):
    node_y = j * dy
    if j % 2 == 0:
        for i in range(-nx + 1, nx):
            node_x = i * dx
            # Exclude points in crack region
            if not (node_x <= crack_x and (crack_y-2*d_equil_y < node_y < crack_y+2*d_equil_y)):
                atom_id = atom_id + 1
                fw1.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(atom_id), str(node_x), str(node_y), str(dx)))
                fw3.write("%s 1 %s %s 0.000000e+00 0 0 0\n" % (str(atom_id), str(node_x), str(node_y)))
                if node_y == crack_y+2*d_equil_y and node_x <= crack_x:
                    id_atom_crack_upper.append(atom_id)
                elif node_y == crack_y-2*d_equil_y and node_x <= crack_x:
                    id_atom_crack_lower.append(atom_id)
                elif node_y == crack_y and node_x == crack_x+d_equil:
                    id_atom_crack_tip.append(atom_id)
                else:
                    id_atom_point_inside.append(atom_id)

    else:
        for i in range(-nx + 1, nx - 1):
            node_x = dx / 2 + i * dx
            if not (node_x <= crack_x and (crack_y-2*d_equil_y < node_y < crack_y+2*d_equil_y)):
                atom_id = atom_id + 1
                fw1.write("Point(%s) = {%s, %s, 0.0, %s};\n" % (str(atom_id), str(node_x), str(node_y), str(dx)))
                fw3.write("%s 1 %s %s 0.000000e+00 0 0 0\n" % (str(atom_id), str(node_x), str(node_y)))
                if (crack_x < node_x < crack_x+d_equil) and node_y == crack_y+d_equil_y:
                    id_atom_crack_tip_upper.append(atom_id)
                elif (crack_x < node_x < crack_x+d_equil) and node_y == crack_y-d_equil_y:
                    id_atom_crack_tip_lower.append(atom_id)
                else:
                    id_atom_point_inside.append(atom_id)

# ------------------------ generate mesh for detection band ------------------------------

for i in range (0, len(id_1234_inner)-1):
    line_id = id_1234_inner[i]
    fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_1234_inner[i]), str(id_1234_inner[i+1])))

line_id = line_id + 1
fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_1234_inner[-1]), str(id_5_crack_upper[0])))

for i in range(-1, len(id_atom_crack_upper)-1):
    line_id = line_id + 1
    if i == -1:
        fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_5_crack_upper[0]), str(id_atom_crack_upper[0])))
    else:
        fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_atom_crack_upper[i]), str(id_atom_crack_upper[i+1])))

line_id = line_id + 1
fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_atom_crack_upper[-1]), str(id_atom_crack_tip_upper[0])))
line_id = line_id + 1
fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_atom_crack_tip_upper[0]), str(id_atom_crack_tip[0])))
line_id = line_id + 1
fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_atom_crack_tip[0]), str(id_atom_crack_tip_lower[0])))

for i in range(0, len(id_atom_crack_lower)):
    line_id = line_id + 1
    if i == 0:
        fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_atom_crack_tip_lower[0]), str(id_atom_crack_lower[-1])))
    else:
        fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_atom_crack_lower[-i]), str(id_atom_crack_lower[-i-1])))

line_id = line_id + 1
fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_atom_crack_lower[0]), str(id_11_crack_lower[0])))

for i in range(-1, len(id_12_inner)-1):
    line_id = line_id+1
    if i == -1:
        fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_11_crack_lower[0]), str(id_12_inner[0])))
        phy_line_inner2_min = line_id
    else:
        fw1.write("Line(%s) = {%s, %s};\n" % (str(line_id), str(id_12_inner[i]), str(id_12_inner[i+1])))

last_line_id = line_id+1
fw1.write("Line(%s) = {%s, %s};\n" % (str(last_line_id), str(id_12_inner[-1]), str(id_1234_inner[0])))

fw1.write("Curve Loop(1) = {1:%s};\n" % str(last_line_id))
fw1.write("Plane Surface(1) = {1};\n")

for i in range(0, len(id_atom_point_inside)):
    fw1.write("Point {%s} In Surface {1};\n" % (str(id_atom_point_inside[i])))

# ------------------------------------------------------------------------------

fw.close()
fw1.close()
fw3.close()



