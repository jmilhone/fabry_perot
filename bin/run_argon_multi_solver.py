from __future__ import print_function, division
import numpy as np
import pymultinest
from mpi4py import MPI
import argparse
import os.path as path
from distutils.util import strtobool
import time
from fabry.tools import file_io


def parse_csv_config(filename):
    """Parses csv config file since pandas can't be installed on Wisc HPC

    :param str filename: path to config file to be parsed
    :returns: a dictionary containing the configuration information, everything is a string
    :rtype: dict
    """
    info_dict = dict()
    with open(filename, 'r') as infile:
        # remove newline, and split on delimiter ','
        header = infile.readline().rstrip().split(',')

        # build dict with header names as the keys
        for name in header:
            info_dict[name] = list()

        for line in infile:
            # remove newline, and split on delimiter ','
            line_data = line.rstrip().split(",")

            # iterate over header values and data and add to info_dict[key]
            for key, value in zip(header, line_data):
                info_dict[key].append(value)
    return info_dict

def convert_items_to_type(input_dict, keys_to_convert, types):
    """Converts values in input_dict[key] for each key in keys_to_convert to the type in types

    :param dict input_dict: configuration dictionary to modify
    :param list keys_to_convert: keys in input_dict to convert value types
    :param list types: list of type callable functions to convert each key value with
    """
    for key, f in zip(keys_to_convert, types):
        value = input_dict[key]
        value = [f(x) for x in value]
        input_dict[key] = value

def retrieve_ringsum_data(filename, image_index):
    """Reads the ringsum data at image_index in the hdf5 filename

    :param filename: path to file to read in ringsum data
    :param int image_index: index from the different time points from a plasma shot to read
    :returns: Tuple with radius, signal, and signal std arrays
    :rtype: Tuple
    """
    data = file_io.h5_2_dict(filename)

    r = data['r'][image_index]
    sig = data['sig'][image_index]
    sd = data['sd'][image_index]

    return r, sig, sd

def chop_down_to_fit_range(r, signal, std, fit_range):
    """Reduces the portion of the ringsum to fit_range

    :param r: ringsum radius array
    :param signal: ringsum singal array
    :param std: ringsum standard error array
    :param fit_range: lower and upper bound in the radius array to chop down to
    :returns: Tuple with chopped down radius, signal, and signal std arrays
    """
    i0 = np.abs(r-fit_range[0]).argmin()
    i1 = np.abs(r-fit_range[1]).argmin()
    sl = slice(i0, i1+1) 

    return r[sl], signal[sl], std[sl]



def read_calibration(calibration_folder):
    """
    """
    posterior = np.loadtxt(path.join(calibration_folder, 'full_post_equal_weights.dat'), ndmin=2)
    # L, d, F are indices 0, 1, and 2
    return posterior[:, 0:3]


def verify_restart():
    """
    """
    #a = raw_input("Are you sure you want to restart? ")
    a = input("Are you sure you want to restart? ")
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
    """
    """
    clean_path = path.expanduser(filepath)
    clean_path = path.abspath(clean_path)
    return clean_path

def main():
    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()

    if rank == 0:
        # single process here
        parser = argparse.ArgumentParser()
        parser.add_argument('config', type=str,
                help='Filepath of configuration file')
        parser.add_argument('cal_folder', metavar='calibration-folder', type=str,
                help='Folder containing the fabry-perot calibration')
        parser.add_argument('out_folder', metavar='output-folder',
                type=str, help='Folder to place output files in')
        parser.add_argument('--restart', action='store_true', default=False,
                            help="Set to True if you want MultiNest to start over instead of resuming")
        parser.add_argument('--image_index', type=int, default=5,
                help="Image index to solve for")

        args = parser.parse_args()

        restart = args.restart
        if restart:
            restart = verify_restart()

        output_folder = clean_filepath(args.out_folder)
        config_name = clean_filepath(args.config)
        cal_folder = clean_filepath(args.cal_folder)

        calib_posterior = read_calibration(cal_folder)

        info_dict = parse_csv_config(config_name) 
        convert_items_to_type(info_dict, ['impact_factor', 'gain', 'int_time', 'shot_num'], types=[float, float, float, int])

        image_index = args.image_index
        filepaths = info_dict['filepath']
        filenames = info_dict['filename']
        filenames = [path.join(fpath, fname) for (fpath, fname) in zip(filepaths, filenames)]

        impact_factors = info_dict['impact_factor']

        r = list()
        sig = list()
        sd = list()
        # grab the ringsum data at index=image_index
        for fname in filenames:
            data = retrieve_ringsum_data(fname, image_index)
            r.append(data[0])
            sig.append(data[1])
            sd.append(data[2])

        n_images = len(r)

        fit_range = (193, 230)
        # Only look at the data if r is in interval [fit_range[0], fit_range[1]]
        for idx in range(n_images):
            r[idx], sig[idx], sd[idx] = chop_down_to_fit_range(r[idx], sig[idx], sd[idx], fit_range)

        # Adjust signal/std amplitudes based on camera gain and integration time 
        for idx, (signal, std) in enumerate(zip(sig, sd)):
            #signal *= info_df.loc[idx, 'gain'] / info_df.loc[idx, 'int_time']
            #std *= info_df.loc[idx, 'gain'] / info_df.loc[idx, 'int_time']
            # moving toward a dict instead of a dataframe
            signal *= info_dict['gain'][idx] / info_dict['int_time'][idx]
            std *= info_dict['gain'][idx] / info_dict['int_time'][idx]

        solver_input = {'output_folder': output_folder,
                        'calib_posterior': calib_posterior,
                        'restart': restart,
                        'config_name': config_name,
                        # 'image_index': args.image_index,
                        'r': r,
                        'sig': sig,
                        'sd': sd,
                        'impact_factors': impact_factors
                       }
    else:
        solver_input = None

    # important line here, this is the brodcasting to the other processes
    solver_input = Comm.bcast(solver_input, root=0)

    if solver_input is not None:
        from fabry.plasma import argon_chord_solver as acs

        # Expand solver_in for readability
        output_folder = solver_input['output_folder']
        calib_posterior = solver_input['calib_posterior']
        restart = solver_input['restart']
        config_name = solver_input['config_name']
        # image_index = solver_input['image_index']

        image_data = (solver_input['r'], solver_input['sig'], solver_input['sd'], solver_input['impact_factors'])
        # run the solver here
        start_time = time.time()
        print("solving...")
        acs.argon_multi_image_solver_fixed_Lnu(output_folder, calib_posterior, config_name, image_data,
                                               resume=not restart, test_plot=False)
        end_time = time.time()

    if rank == 0:
        print("Total Time Elapsed: {} minutes".format((end_time-start_time)/60.0))
        # can also run the posterior plotting here...
        from fabry.plasma import check_argon_chord_solver as cacs
        cacs.check_const_Lnu_solver(output_folder, calib_posterior, image_data)
if __name__ == "__main__":
    main()
