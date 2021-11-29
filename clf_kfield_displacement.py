# Coupled LAMMPS and FEniCS
# To calculate the displacement bc to be applied based on the prescribed stress intensity factor
# Written by Wenjia G. on 08/06/2018

import math

# =============== Displacement bc functions ============== #
def kfield_disp(x, y, crack_x, crack_y, K_I, E, nu, kappa):
    pi = math.pi

    r = math.sqrt((x-crack_x)**2+(y-crack_y)**2)
    theta = math.atan2(y-crack_y,x-crack_x)

    ux = K_I/(2*E)*math.sqrt(r/(2*math.pi))*((1+nu)*((2*kappa-1)*math.cos(theta/2)-math.cos(3*theta/2)))
    uy = K_I/(2*E)*math.sqrt(r/(2*math.pi))*((1+nu)*((2*kappa+1)*math.sin(theta/2)-math.sin(3*theta/2)))

    return ux, uy



