from __future__ import division, print_function
import sys
from future.builtins import input
from distutils.util import strtobool
import argparse
import os.path as path
from mpi4py import MPI
sys.path.append("../")
from fabry.tools import file_io
import numpy as np

def verify_restart():
    a = raw_input("Are you sure you want to restart? ")
    try:
        a = bool(strtobool(a))
    except ValueError:
        print("invalid input exiting...")
        sys.exit(1)

    if a:
        print("ok moving on with restart")
        restart = True
    else:
        print("Ok, overwriting restart to False")
        restart = False

    return restart


if __name__ == "__main__":
    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    if rank == 0:
        parser = argparse.ArgumentParser(description='Performs a Ti and V solver in Ar.')

        parser.add_argument('folder', type=str,
                help='folder of finesse image to use, this is also where multinest data is saved')
        parser.add_argument('ld_folder', type=str, help='folder containing Ld calibration multinest output')
        parser.add_argument('finesse_folder', type=str,
                help='folder containing finesse calibration multinest output')
        parser.add_argument('filter_type', type=str,
                            help='Name of filter (Ar, He, ar, he, ...)')
        parser.add_argument('--restart', action='store_true', default=False,
                            help="Set to True if you want MultiNest to start over instead of resuming")
        args = parser.parse_args()

        if args.filter_type.lower() in ['ar', 'argon', '488']:
            filter_type = 'argon'
        elif args.filter_type.lower() in ['he', 'helium', '468.6']:
            filter_type = 'helium'
        else:
            raise NotImplementedError('Filter {0:s} not recognized.'.format(args.filter_type))

        folder = path.abspath(path.expanduser(args.folder))
        ld_folder = path.abspath(path.expanduser(args.ld_folder))
        finesse_folder = path.abspath(path.expanduser(args.finesse_folder))

        # L_posterior, d_posterior = file_io.read_Ld_results(ld_folder)
        #F_posterior, _, _, _, _ = file_io.read_finesse_results(finesse_folder)

        Ldpost = np.loadtxt(path.join(path.abspath(args.ld_folder), 'full_post_equal_weights.dat'), ndmin=2)
        L_posterior = Ldpost[:, 0]
        d_posterior = Ldpost[:, 1]
        F_posterior = Ldpost[:, 2]
        restart = args.restart
        if restart:
            restart = verify_restart()
        solver_in = {'folder': folder,
                     'Lpost': L_posterior,
                     'dpost': d_posterior,
                     'Fpost': F_posterior,
                     'restart': restart,
                     'filter': filter_type,
                     }

    else:
        solver_in = None

    solver_in = Comm.bcast(solver_in, root=0)

    if solver_in is not None:
        # here is where we will call the solver
        if solver_in['filter'] == 'argon':
            # call argon solver of choice
            import fabry.plasma.argon_plasma_solver as solver
            solver_args = (solver_in['folder'], solver_in['Fpost'], solver_in['Lpost'], solver_in['dpost'])

            #solver.no_vel_solver(*solver_args, resume=not solver_in['restart'], test_plot=True)
            #solver.const_vel_solver(*solver_args, resume=not solver_in['restart'], test_plot=False)
            #solver.profile_vel_solver(*solver_args, resume=not solver_in['restart'], test_plot=False)

        elif solver_in['filter'] == 'helium':
            raise NotImplementedError("Helium is not implemented yet")
        else:
            raise ValueError("How did we get here?")


    if rank == 0:
        # here is where we will check the solver
        if solver_in['filter'] == 'argon':
            # call argon solver of choice
            import fabry.plasma.check_argon_solver as checker

            solver_args = (solver_in['folder'], solver_in['Fpost'], solver_in['Lpost'], solver_in['dpost'])

            #checker.check_no_solver(*solver_args)
            #checker.check_constV_solver(*solver_args)
            checker.check_profile_solver(*solver_args)

        elif solver_in['filter'] == 'helium':
            raise NotImplementedError("Helium is not implemented yet")
        else:
            raise ValueError("How did we get here?")




