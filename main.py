# Coupled LAMMPS and FEniCS
# The main loop - with multiprocessing package imported
# Written by Wenjia G. om 10/09/2018


# =====================================================

def lammps_run(n_cycle, icycle_p, maxmin_p, convergence_p, updated_atom_p, directory, atomdata, natom, xlo, xhi, ylo, yhi):

    from lammps import lammps
    import os
    import glob
    import numpy as np
    import ctypes
    import shutil

    lammps_step = 0

    while True:

        cycle = icycle_p.recv()
        maxormin = maxmin_p.recv()
        ifconvergence = convergence_p.recv()
        atom_data_new = updated_atom_p.recv()

        if cycle == 0:

            lmp = lammps()

            lmp.file("clf_lammps.inp")

            shutil.move('dump.min_0.cfg', directory+'/'+"atoms.0.cfg")
            for file in glob.glob("dump.*.cfg"):
                os.remove(file)

        else:

            if ifconvergence == 1:
                lammps_step = lammps_step + 1

            c_atom_data_new = (ctypes.c_double * len(atom_data_new))(*atom_data_new)
            lmp.scatter_atoms("x", 1, 3, c_atom_data_new)

            if (maxormin == 1 or maxormin == 2 or lammps_step % 5 == 0) and ifconvergence == 1:
                lmp.command("dump       2 all cfg 10 dump.min_*.cfg mass type xsu ysu zsu vx vy vz")
            lmp.command("minimize       0.0 1.0e-3 10 1000")
            if (maxormin == 1 or maxormin == 2 or lammps_step % 5 == 0) and ifconvergence == 1:
                lmp.command("undump     2")

            for file in glob.glob("dump.*.cfg"):
                shutil.move(file, directory+'/'+"atoms.%s.%s.cfg" % (str(cycle), str(lammps_step)))

        data_0 = lmp.gather_atoms("x", 1, 3)
        data = np.ctypeslib.as_array(data_0) # !! Convert data type for later computation

        atomdata.send(data)
        natom.send(lmp.extract_global("natoms", 0))
        xlo.send(lmp.extract_global("boxxlo", 1))
        xhi.send(lmp.extract_global("boxxhi", 1))
        ylo.send(lmp.extract_global("boxylo", 1))
        yhi.send(lmp.extract_global("boxyhi", 1))

        if cycle == n_cycle and maxormin == 2 and ifconvergence == 1:
            break


# =====================================================


if __name__ == "__main__":

    import multiprocessing as mp
    from clf_key_word_search import *
    from clf_atoms_for_lammps import *
    from clf_fenics_solve import fenics_initial_setup, fenics_displacement
    from clf_crack_tip_function import crack_tip_finding
    import ctypes
    import time
    import numpy as np
    import sys

    sys.setrecursionlimit(200000)

    t_start = time.time()
    t_lammps = 0
    t_lammps_data = 0
    t_fenics = 0

    mp.set_start_method('spawn')

    # =================== Parameter setup ===============

    t_setup_start = time.time()

    input = sys.stdin.readlines()

    trial_id = str(key_search('trial', input)) # Run id
    delta_K_I = float(key_search('delta_K_I', input))  # max stress intensity factor - min stress intensity factor
    r_value = float(key_search('R_value', input)) # min stress intensity factor / max stress intensity factor
    n_cycle = int(key_search('cycles', input)) # Number of cycles
    increment = float(key_search('increment', input)) # loading/unloading increment

    # Read in elastic constants
    E = float(key_search('modulus', input)) # 3D Young's modulus
    nu = float(key_search('ratio', input))  # 3D Poisson's ratio

    mu = E/(2*(1+nu))  # shear modulus
    lmbda = E*nu/((1+nu)*(1-2*nu)) # Lame's first parameter
    kappa = 3 - 4 * nu

    # Read in the prepared data file names
    mesh_filename = str(key_search('mesh', input))
    atom_filename = str(key_search('atom_data', input))

    # Read in initial crack tip position
    x_crack_tip = float(key_search('position x', input))
    y_crack_tip = float(key_search('position y', input))

    # Create a directory for outputs
    directory = "{}-outputs".format(trial_id)
    if not os.path.exists(directory):
        os.mkdir(directory)

    K_I = 0.0
    K_max = delta_K_I/(1-r_value)
    K_min = r_value*K_max


    # =================== Initial setup ===============

    # creating condition object
    condition = mp.Condition()

    # creating variables
    icycle_p, icycle_c = mp.Pipe()
    maxmin_p, maxmin_c = mp.Pipe()
    atomdata_p, atomdata_c = mp.Pipe()
    natom_p, natom_c = mp.Pipe()
    xlo_p, xlo_c = mp.Pipe()
    xhi_p, xhi_c = mp.Pipe()
    ylo_p, ylo_c = mp.Pipe()
    yhi_p, yhi_c = mp.Pipe()
    convergence_p, convergence_c = mp.Pipe()
    updated_atom_p, updated_atom_c = mp.Pipe()

    # creating new process
    p_lammps = mp.Process(target=lammps_run, args=(n_cycle, icycle_p, maxmin_p, convergence_p, updated_atom_p,
                                                   directory, atomdata_c, natom_c, xlo_c, xhi_c, ylo_c, yhi_c))

    # starting process
    p_lammps.start()

    # =~=~=~=~=~=~=~= FEniCS and LAMMPS setup =~=~=~=~=~=~=
    # ------------- Initial FEniCS setup -------------
    mesh, fd, V, d, L, coordinates, v2d, global_K = fenics_initial_setup(mesh_filename, K_I, mu, nu, lmbda, E, kappa, x_crack_tip, y_crack_tip)

    # ---------- Initial LAMMPS setup and first run ----------
    # prepare atom data
    interface_coord, n_atom, pad_atom_coord = lammps_intial_atom_data(atom_filename, x_crack_tip, y_crack_tip, K_I, E, nu, kappa)

    # -----------------------------------------------
    # !!! Send variables and invoke LAMMPS
    with condition:
        icycle_c.send(0)
        maxmin_c.send(0)
        convergence_c.send(0)
        updated_atom_c.send(0)

    # -----------------------------------------------

    atom_data = atomdata_p.recv()
    n_atom = natom_p.recv()
    x_lo = xlo_p.recv()
    x_hi = xhi_p.recv()
    y_lo = ylo_p.recv()
    y_hi = yhi_p.recv()

    last_u_vertex = 0.0  # used for calculating convergence

    t_setup_end = time.time()
    t_setup = t_setup_end - t_setup_start


    # ========================= THE MAIN LOOP =========================

    for icycle in range(1, n_cycle+1):

        #================================ New Cycle Starts! ================================
        # ------------------------ LOADING ------------------------
        while 1:  # loading

            # Crack tip position updated ever loading step
            x_crack_tip, y_crack_tip = crack_tip_finding(n_atom, atom_data, x_lo, x_hi, y_lo, y_hi)

            with open("crack_tip_position.txt", "a") as myfile:
                myfile.write("------ Cycle # %s : K_I = %s ------ " % (str(icycle), str(K_I)))
                myfile.write('\n')
                myfile.write("%s %s" % (str(x_crack_tip), str(y_crack_tip)))
                myfile.write('\n')

            maxmin = 0
            ifconvergence = 0

            if K_I >= K_max: 
                maxmin = 1 

            istep = 0  # Loop until result converges

            while 1:
                istep = istep + 1

                # ------------- Compute displacements in FEniCS -------------
                t_fenics_start=time.time()
                last_u_vertex, u_node_max = fenics_displacement(trial_id, icycle, istep, maxmin, atom_data, K_I, mu, nu, lmbda, E,
                                                                kappa, x_crack_tip, y_crack_tip, mesh, fd, V, d, L, v2d,
                                                                global_K, interface_coord, last_u_vertex)
                t_fenics_end=time.time()
                t_fenics = t_fenics+(t_fenics_end-t_fenics_start)

                # ------------- Update pad atom positions and run LAMMPS -------------
                t_lammps_start = time.time()

                # Update pad atom data
                updated_atom = lammps_pad_update(pad_atom_coord, atom_data, last_u_vertex)

                if u_node_max < 1e-4:  # Convergence
                    ifconvergence = 1

                # -----------------------------------------------
                # !!! Send variables and invoke LAMMPS
                with condition:
                    icycle_c.send(icycle)
                    maxmin_c.send(maxmin)
                    convergence_c.send(ifconvergence)
                    updated_atom_c.send(updated_atom)

                # -----------------------------------------------

                atom_data = atomdata_p.recv()
                n_atom = natom_p.recv()
                x_lo = xlo_p.recv()
                x_hi = xhi_p.recv()
                y_lo = ylo_p.recv()
                y_hi = yhi_p.recv()

                t_lammps_end = time.time()
                t_lammps = t_lammps+(t_lammps_end-t_lammps_start)

                if ifconvergence == 1:
                    break

            if maxmin == 1:
                K_I = K_I - increment
                break

            K_I = K_I+increment


        # ------------------------ UNLOADING ------------------------
        while 1:  # unloading

            x_crack_tip, y_crack_tip = crack_tip_finding(n_atom, atom_data, x_lo, x_hi, y_lo, y_hi)

            with open("crack_tip_position.txt", "a") as myfile:
                myfile.write("------ Cycle # %s : K_I = %s ------ " % (str(icycle), str(K_I)))
                myfile.write('\n')
                myfile.write("%s %s" % (str(x_crack_tip), str(y_crack_tip)))
                myfile.write('\n')

            maxmin = 0
            ifconvergence = 0

            if K_I <= K_min:  
                maxmin = 2 

            istep = 0  # Loop until result converges

            while 1:
                istep = istep + 1

                # ------------- Compute displacements in FEniCS -------------
                t_fenics_start = time.time()
                last_u_vertex, u_node_max = fenics_displacement(trial_id, icycle, istep, maxmin, atom_data, K_I, mu, nu, lmbda, E,
                                                                kappa, x_crack_tip, y_crack_tip, mesh, fd, V, d, L, v2d, global_K,
                                                                interface_coord, last_u_vertex)
                t_fenics_end = time.time()
                t_fenics = t_fenics + (t_fenics_end - t_fenics_start)

                # ------------- Update pad atom positions and run LAMMPS -------------
                t_lammps_start = time.time()

                # Update pad atom data
                updated_atom = lammps_pad_update(pad_atom_coord, atom_data, last_u_vertex)

                if u_node_max < 1e-4:  # Convergence
                    ifconvergence = 1

                # -----------------------------------------------
                # !!! Send variables and invoke LAMMPS
                with condition:
                    icycle_c.send(icycle)
                    maxmin_c.send(maxmin)
                    convergence_c.send(ifconvergence)
                    updated_atom_c.send(updated_atom)

                # -----------------------------------------------

                atom_data = atomdata_p.recv()
                n_atom = natom_p.recv()
                x_lo = xlo_p.recv()
                x_hi = xhi_p.recv()
                y_lo = ylo_p.recv()
                y_hi = yhi_p.recv()

                t_lammps_end = time.time()
                t_lammps = t_lammps + (t_lammps_end - t_lammps_start)

                # if ifconvergence == 1 and maxmin == 1:
                if ifconvergence == 1:
                    break

            if maxmin == 2:
                K_I = K_I + increment
                break

            K_I = K_I - increment

       
    p_lammps.join()

    with open("DONE", "w") as myfile:
        myfile.write("ALL DONE!!!")

    t_end = time.time()
    t_total = t_end - t_start

    with open("time", "w") as myfile:
        myfile.write("t_total = %s \n" % str(t_total))
        myfile.write("t_fenics = %s \n" % str(t_fenics))
        myfile.write("t_lammps = %s \n" % str(t_lammps))
        myfile.write("t_lammps_data = %s \n" % str(t_lammps_data))
        myfile.write("t_setup = %s \n" % str(t_setup))






