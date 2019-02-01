from __future__ import division, print_function
#import sys
import os.path as path
from distutils.util import strtobool
from mpi4py import MPI
#sys.path.append("../")
from fabry.tools import file_io
import numpy as np
import argparse
import ConfigParser
import time

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


def clean_filepath(filepath):
    clean_path = path.expanduser(filepath)
    clean_path = path.abspath(clean_path)
    return clean_path


def read_calibration(calibration_folder):
    posterior = np.loadtxt(path.join(calibration_folder, 'full_post_equal_weights.dat'), ndmin=2)
    L_posterior = posterior[:, 0]
    d_posterior = posterior[:, 1]
    F_posterior = posterior[:, 2]
    return L_posterior, d_posterior, F_posterior


def parse_config(config_filename):
    config = ConfigParser.ConfigParser()
    config.read(clean_filepath(config_filename))
    sections = config.sections()

    locs = []
    folders = []

    for section in sections:
        locs.append(float(section))
        folders.append(clean_filepath(config.get(section, 'path'))) 
    return locs, folders

if __name__ == "__main__":
    start_time = time.time()
    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    if rank == 0:
        parser = argparse.ArgumentParser(description='Performs a Ti and PCX velocity profile solver in Ar using multiple images.')
        parser.add_argument('config', type=str,
                help='Filepath of configuration file')
        parser.add_argument('cal_folder', metavar='calibration-folder', type=str,
                help='Folder containing the fabry-perot calibration')
        parser.add_argument('out_folder', metavar='output-folder',
                type=str, help='Folder to place output files in')
        parser.add_argument('filter_type', metavar='filter-type',type=str,
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

        restart = args.restart
        if restart:
            restart = verify_restart()

        output_folder = clean_filepath(args.out_folder)

        config_filename = clean_filepath(args.config)
        chord_locs, folders = parse_config(config_filename)

        for loc, folder in zip(chord_locs, folders):
            print(loc, folder)

        cal_folder = clean_filepath(args.cal_folder)
        L_post, d_post, F_post = read_calibration(cal_folder)

        solver_in = {'Lpost': L_post,
                     'dpost': d_post,
                     'Fpost': F_post,
                     'chord_locs': chord_locs,
                     'folders': folders,
                     'output_folder':output_folder,
                     'restart': restart,
                     'filter': filter_type,
                     }
    else:
        solver_in = None

    solver_in = Comm.bcast(solver_in, root=0)

    if solver_in is not None:
        # run everything
        import fabry.plasma.argon_plasma_solver as solver
        Lpost = solver_in['Lpost']
        dpost = solver_in['dpost']
        Fpost = solver_in['Fpost']
        chord_locs = solver_in['chord_locs']
        folders = solver_in['folders']
        output_folder = solver_in['output_folder']
        restart = solver_in['restart']

        if solver_in['filter'] == 'argon':
            solver.multi_image_solver(output_folder, chord_locs, folders, Lpost, dpost, Fpost, test_plot=False, resume=not restart)
        elif solver_in['filter'] == 'helium':
            raise NotImplementedError("Helium is not implemented yet")
        else:
            raise ValueError("How did we get here?")

    if rank == 0:
        end_time = time.time()
        print("Total Time Elasped: {} minutes".format((end_time-start_time)/60.0))
    #     # run plasma solver check here
    #     if solver_in['filter'] == 'argon':
    #         import fabry.plasma.check_argon_solver as checker
    #         checker.check_multi_solver(output_folder, chord_locs, folders, Lpost, dpost, Fpost)
    #     elif solver_in['filter'] == 'helium':
    #         raise NotImplementedError("Helium is not implemented yet")
    #     else:
    #         raise ValueError("How did we get here?")




