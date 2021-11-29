# Coupled LAMMPS and FEniCS
# To find the position of the crack tip
# Adapted from Rick's work
# Written by Wenjia G. on 08/21/2018

import math
import numpy as np

# ===============================================================
def floodfill(matrix, x, y, y_stop1, y_stop2):

    if matrix[x][y] == 0:
        matrix[x][y] = 2

        # recursively invoke flood fill on all surrounding cells:
        if y_stop1 <= y <= y_stop2:
            if x > 1:
                floodfill(matrix, x-1, y, y_stop1, y_stop2)
                floodfill(matrix, x-1, y-1, y_stop1, y_stop2)
                floodfill(matrix, x-1, y+1, y_stop1, y_stop2)
            if x < len(matrix)-1:
                floodfill(matrix, x+1, y, y_stop1, y_stop2)
                floodfill(matrix, x+1, y-1, y_stop1, y_stop2)
                floodfill(matrix, x+1, y+1, y_stop1, y_stop2)
# ===============================================================

def crack_tip_finding(n_atom, atoms_data, xlo, xhi, ylo, yhi):

    # ------- input parameter -------
    grid_dx = 1.2
    grid_dy = 1.2

    grid_xmin = xlo
    grid_xmax = xhi
    grid_ymin = ylo
    grid_ymax = yhi

    nx_grid = int((grid_xmax-grid_xmin)/grid_dx)+1
    ny_grid = int((grid_ymax-grid_ymin)/grid_dy)+1

    # Create matrix to record if there is any atom in one grid
    CV_cell_free = np.zeros((nx_grid,ny_grid))
    n_CV_cells = nx_grid*ny_grid

    # ------- Assign atoms to different grids -------
    for i in range(n_atom):
        atom_x = atoms_data[3*i]
        atom_y = atoms_data[3*i + 1]

        xsite_d = (atom_x-grid_xmin)/(nx_grid*grid_dx)
        xsite_i = int(math.ceil(xsite_d*nx_grid))
        ysite_d = (atom_y-grid_ymin)/(ny_grid*grid_dy)
        ysite_i = int(math.ceil(ysite_d*ny_grid))

        CV_cell_free[xsite_i-1,ysite_i-1] = 1
        #print(CV_cell_free)

    # ------- Use flood fill algorithm to color the crack -------
    start_i = 1
    start_j = int(ny_grid/2)
    stop_j1 = int(ny_grid/20)
    stop_j2 = ny_grid - int(ny_grid/20)
    floodfill(CV_cell_free, start_i, start_j, stop_j1, stop_j2)

    if np.count_nonzero(CV_cell_free == 2) <= 3:
        for nx in range(1, 11):
            start_i = 1 + nx
            floodfill(CV_cell_free, start_i, start_j, stop_j1, stop_j2)
            if np.count_nonzero(CV_cell_free == 2) > 3:
                break
            else:
                for ny in range(1, 3):
                    for ny_2 in range(2):
                        start_i = 1
                        start_j = int(ny_grid / 2) + (-1) ** ny_2 * ny
                        floodfill(CV_cell_free, start_i, start_j, stop_j1, stop_j2)
                        if np.count_nonzero(CV_cell_free == 2) > 3:
                            break
                        else:
                            for nx in range(1, 11):
                                start_i = 1 + nx
                                floodfill(CV_cell_free, start_i, start_j, stop_j1, stop_j2)
                                if np.count_nonzero(CV_cell_free == 2) > 3:
                                    break

    #print(CV_cell_free)

    # ------- Find the furthest cell that has the value of 2  -------
    max_i = 0
    max_j = 0

    for i in range(1, nx_grid):
        for j in range(int(ny_grid/2)-int(ny_grid/10), int(ny_grid/2)+int(ny_grid/10)):
            if CV_cell_free[i,j] == 2:
                max_i = i
                max_j = j
    #print(max_i, max_j)

    # ------- Correct the bias in y direction -------
    j = max_j
    j_count = 1
    while CV_cell_free[max_i,j-1] == 2:
        j = j-1
        j_count = j_count+1
    y_correction = j_count*grid_dy/2
    #print(j_count)

    x_cracktip_CV = grid_xmin+(max_i+1)*grid_dx
    y_cracktip_CV = grid_ymin+(max_j+1)*grid_dy-y_correction
    #print(x_cracktip_CV, y_cracktip_CV)

    return x_cracktip_CV, y_cracktip_CV


