# Coupled LAMMPS and FEniCS
# To FEniCS part
# Written by Mingjie Z. 
# Modified by Wenjia G. on 09/25/2018

from __future__ import print_function
from fenics import *
import numpy as np
import math
import sys
from clf_key_word_search import *


pi = math.pi

# ============== Define strain and stress ==============
# plane strain

def epsilon(u):
    return sym(nabla_grad(u))
    # return 0.5*(nabla_grad(u) + nabla_grad(u).T)

def sigma(u,mu,lmbda,d):
    return lmbda*tr(epsilon(u))*Identity(d) + 2.0*mu*epsilon(u)
    # return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)

# ============================ Initial FEniCS setup ============================
def fenics_initial_setup(mesh_filename, K_I, mu, nu, lmbda, E, kappa, crack_x, crack_y):

    # ------------ Initial FEniCS setup ------------

    mu = Constant(mu)
    nu = Constant(nu)
    lmbda = Constant(lmbda)
    E = Constant(E)
    kappa = Constant(kappa)

    # ----------- Create mesh and function space -----------

    mesh_names = mesh_filename.split('.')
    meshcd_filename = mesh_names[0]+'_physical_region.xml'
    meshfd_filename = mesh_names[0]+'_facet_region.xml'

    mesh = Mesh(mesh_filename)
    cd = MeshFunction('size_t', mesh, meshcd_filename);
    fd = MeshFunction('size_t', mesh, meshfd_filename);

    V = VectorFunctionSpace(mesh, 'P', 1)

    # ---------------  Boundary Condition ---------------
    u_D = Expression((
        'K_I/(2*E)*sqrt((sqrt(pow(x[0]-crack_x,2)+pow(x[1]-crack_y,2)))/(2*pi))*((1+nu)*((2*kappa-1)*cos((atan2(x[1]-crack_y,x[0]-crack_x))/2)-cos(3*(atan2(x[1]-crack_y,x[0]-crack_x))/2)))',
        'K_I/(2*E)*sqrt((sqrt(pow(x[0]-crack_x,2)+pow(x[1]-crack_y,2)))/(2*pi))*((1+nu)*((2*kappa+1)*sin((atan2(x[1]-crack_y,x[0]-crack_x))/2)-sin(3*(atan2(x[1]-crack_y,x[0]-crack_x))/2)))'),
        degree=1, K_I=K_I, mu=mu, nu=nu, crack_x=crack_x, crack_y=crack_y, pi=pi, E=E, kappa=kappa)

    bc_outer = DirichletBC(V, u_D, fd, 22)
    bc_interface = DirichletBC(V, u_D, fd, 11)

    bcs = [bc_outer, bc_interface]

    # -------------- Define variational problem ------------
    u = TrialFunction(V)
    d = u.geometric_dimension()  # space dimension = 2
    v = TestFunction(V)
    f = Constant((0, 0))
    T = Constant((0, 0))
    a = inner(sigma(u,mu,lmbda,d), epsilon(v))*dx
    L = dot(f, v)*dx + dot(T, v)*ds

    # --------------- Explicitly Assemble K -------------
    global_K = assemble(a)
    global_f = assemble(L)

    # ----------------- Apply First B.C. ----------------
    [bc.apply(global_K, global_f) for bc in bcs]

    # --------------- CHECK Vertex vs DOF ---------------

    # Assign global force vector into u to post process
    uu = Function(V)
    for i_i in range(len(global_f.array())):
        uu.vector()[i_i] = global_f.array()[i_i]
    vertex_values_0 = uu.compute_vertex_values()

    # Obtain mesh ordinates and its corresponding vertex values
    coordinates = mesh.coordinates()

    # DOF / Vertex Mapping
    v2d_raw = vertex_to_dof_map(V)
    d2v_raw = dof_to_vertex_map(V)

    v2d = v2d_raw.reshape((-1, mesh.geometry().dim()))
    d2v = d2v_raw[range(0, len(d2v_raw), 2)] / 2

    return mesh, fd, V, d, L, coordinates, v2d, global_K


# ============================ FEniCS computation for displacements ============================
def fenics_displacement(trial_id, icycle, istep, maxmin, atoms_data, K_I, mu, nu, lmbda, E, kappa, crack_x, crack_y,
                        mesh, fd, V, d, L, v2d, global_K, interface_coord, last_u_vertex):

    # ------------ FEniCS compuation for displacements ------------

    mu = Constant(mu)
    nu = Constant(nu)
    lmbda = Constant(lmbda)
    E = Constant(E)
    kappa = Constant(kappa)

    # ------------- Update the outer boundary bc in FEniCS -------------
    u_D = Expression((
        'K_I/(2*E)*sqrt((sqrt(pow(x[0]-crack_x,2)+pow(x[1]-crack_y,2)))/(2*pi))*((1+nu)*((2*kappa-1)*cos((atan2(x[1]-crack_y,x[0]-crack_x))/2)-cos(3*(atan2(x[1]-crack_y,x[0]-crack_x))/2)))',
        'K_I/(2*E)*sqrt((sqrt(pow(x[0]-crack_x,2)+pow(x[1]-crack_y,2)))/(2*pi))*((1+nu)*((2*kappa+1)*sin((atan2(x[1]-crack_y,x[0]-crack_x))/2)-sin(3*(atan2(x[1]-crack_y,x[0]-crack_x))/2)))'),
        degree=1, K_I=K_I, mu=mu, nu=nu, crack_x=crack_x, crack_y=crack_y, pi=pi, E=E, kappa=kappa)

    bc = DirichletBC(V, u_D, fd, 22)
    #bc_interface = DirichletBC(V, u_D, fd, 11)
    #bcs = [bc_outer, bc_interface]

    global_f = assemble(L)

    bc.apply(global_f)

    global_f_1 = global_f.array()

    # ------------- Update the interface nodes bc in FEniCS -------------
    # Create a new force vector
    global_f_updated = []
    for i_f in range(len(global_f_1)):
        global_f_updated.append(global_f_1[i_f])
    global_f_updated = np.asarray(global_f_updated)

    # # Update displacement: Current position minus coordinate
    for i in range(len(interface_coord)):
        nodeA = interface_coord[i][0]
        nodeA_x_position = atoms_data[3 * nodeA - 3]
        nodeA_y_position = atoms_data[3 * nodeA - 2]

        nodeA_x_disp = nodeA_x_position - interface_coord[i][1]
        nodeA_y_disp = nodeA_y_position - interface_coord[i][2]

        # Manually assign Displacements into Force Vector
        # **** Vertex #: starts at 1; DOF #: starts at 0
        dof_nodeA_x = v2d[nodeA - 1][0]
        dof_nodeA_y = v2d[nodeA - 1][1]

        global_f_updated[dof_nodeA_x] = nodeA_x_disp
        global_f_updated[dof_nodeA_y] = nodeA_y_disp


    # ------------------- Compute displacement field -------------------

    u = Function(V)

    f = Function(V)
    f.vector().set_local(global_f_updated)

    solve(global_K, u.vector(), f.vector())

    # Output disp ordered by vertax id
    u_vertex_values = u.compute_vertex_values()

    difference = u_vertex_values - last_u_vertex
    u_node_max = max(difference)

    with open("convergence_study.txt", "a") as myfile:
        if istep == 1:
            myfile.write("------ Cycle # %s : K_I = %s ------ " % (str(icycle), str(K_I)))
            myfile.write('\n')
        myfile.write(str(u_node_max))
        myfile.write('\n')

    last_u_vertex_new = np.copy(u_vertex_values)
                
    return last_u_vertex_new, u_node_max








