from __future__ import print_function, division
import numpy as np
import pymultinest
import os.path
# Go back and look into relative imports! 
from fabry.tools import file_io
import fabry.plasma.plasma as plasma
import fabry.core.models as models
import fabry.distributions.distribution as dist

w0 = 487.98634
mu = 39.948


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


def pcx_vfd_model(r, b, Ti_params, V_params, ne_params, delta_d, amp, L, d, F):
    w_model, radiance_model = plasma.vfd_chord(b, Ti_params, V_params, ne_params,
            rmax=40.0, w0=w0, mu=mu, nr=101, nlambda=2048, test_plot=False)

    model_signal = amp*models.general_model(r, L, d+delta_d, F, w_model, radiance_model)
    return model_signal


# def argon_multi_image_solver_fixed_Lnu(output_folder, calib_posterior, config_name,  image_index=4, resume=True, test_plot=False):
def argon_multi_image_solver_fixed_Lnu(output_folder, calib_posterior, 
        config_name,  image_data, image_index=4, resume=True, test_plot=False):
    """Runs MultiNest solver for multiple chords in a PCX Ar plasma with a fixed Lnu.

    :param str output_folder: output folder to store multinest files
    :param np.ndarray calib_posterior: calibration equally weighted posterior samples
    :param DataFrame info_dict: pandas DataFrame containing information for multichord solver
    :param int image_index: selects which image from orginal stacked tiff to analyze
    :param bool resume: Resume multinest solver if True, default=True
    :param bool test_plot: Runs test plot sampling from prior instead of running solver, default=False
    """

    def prior(cube, ndim, nparams):
        """Transforms prior hypercube to physical units
        """
        cube[0] = dist.UniformPrior(cube[0], 0.025, 3.0)
        cube[1] = dist.UniformPrior(cube[1], -5000.0, 5000.0)

        # I want to use a less broad prior for the amplitudes, but not allow negative amplitude
        for i in range(2, nparams):
            cube[i] = dist.TruncatedNormal(cube[i], 60.0, 30.0, 0.0, 1000.0)
            # cube[i] = dist.UniformPrior(cube[i], 1.0, 100.0)

    def loglikelihood(cube, ndim, nparams):
        """Calculates the log likelihood
        """
        chisq = 0.0

        # Sample from the calibration posterior
        # i_choice = np.random.choice(n_calib_points)
        # L = calib_posterior[i_choice, 0]
        # d = calib_posterior[i_choice, 1]
        # F = calib_posterior[i_choice, 2]

        # For testing purposes, dont' have a calibration to use
        L = 74.6 / 0.00586
        d = 0.8836181
        F = 20.9

        # For each image, calculate the model signal and update the total chisquared
        i = 0
        for loc, rr, ss, std  in zip(impact_factors, r, sig, sd):
            # sample a v_offset to account for the zero velocity calibration drift
            v_offset = np.random.normal(loc=0.0, scale=75.0, size=1)

            wavelength, emission = plasma.calculate_pcx_chord_emission(loc, cube[0],
                    # w0, mu, Lnu, cube[1], rmax=rmax, nr=101, nlambda=2000,
                    w0, mu, Lnu, cube[1]+v_offset, rmax=rmax, nr=101, nlambda=2000,
                    R_outer=r_anode)
            model_signal = models.general_model(rr, L, d, F, wavelength, emission)
            model_signal = cube[i+2] * model_signal

            chisq += np.sum((model_signal-ss)**2 / std**2)

            i+=1

        return -chisq / 2.0

    # info_dict = parse_csv_config(config_name) 
    # convert_items_to_type(info_dict, ['impact_factor', 'gain', 'int_time', 'shot_num'], types=[float, float, float, int])

    # print(info_dict)
    #### Constants
    px_size = 0.00586
    Lnu = 100.0
    rmax = 40.0
    r_anode = 32.0
    fit_range = (180, 230)
    n_calib_points = calib_posterior.shape[0]

    # filepaths = info_dict['filepath']
    # filenames = info_dict['filename']
    # filenames = [os.path.join(fpath, fname) for (fpath, fname) in zip(filepaths, filenames)]

    # impact_factors = info_dict['impact_factor']

    # r = list()
    # sig = list()
    # sd = list()

    # grab the ringsum data at index=image_index
    # for fname in filenames:
    #     data = retrieve_ringsum_data(fname, image_index)
    #     r.append(data[0])
    #     sig.append(data[1])
    #     sd.append(data[2])

    # n_images = len(r)

    # # Only look at the data if r is in interval [fit_range[0], fit_range[1]]
    # for idx in range(n_images):
    #     r[idx], sig[idx], sd[idx] = chop_down_to_fit_range(r[idx], sig[idx], sd[idx], fit_range)

    # # Adjust signal/std amplitudes based on camera gain and integration time 
    # for idx, (signal, std) in enumerate(zip(sig, sd)):
    #     #signal *= info_df.loc[idx, 'gain'] / info_df.loc[idx, 'int_time']
    #     #std *= info_df.loc[idx, 'gain'] / info_df.loc[idx, 'int_time']
    #     # moving toward a dict instead of a dataframe
    #     signal *= info_dict['gain'][idx] / info_dict['int_time'][idx]
    #     std *= info_dict['gain'][idx] / info_dict['int_time'][idx]

    # unpact image_data
    r = image_data[0]
    sig = image_data[1]
    sd = image_data[2]
    impact_factors = image_data[3]
    n_images = len(r)

    # params are Ti, Vmax, +amplitudes for each image
    n_params = 2 + n_images

    if test_plot:
        import matplotlib.pyplot as plt
        # sample from prior and see if things make sense
        nsamples = 500
        L = 74.6 / px_size
        d = 0.8836181
        F = 21.0
        cubes = np.random.rand(nsamples, n_params)

        # grab image data to compare prior samples to
        image_idx = 2
        loc = impact_factors[image_idx]
        rr = r[image_idx]
        ss = sig[image_idx]
        std = sd[image_idx]

        model_signal = np.zeros((nsamples, len(rr)))

        # for idx, cube in enumerate(cubes):
        #     #transform unit hypercube to physical units
        #     prior(cube, None, n_params)
        #     # calculate model
        #     wavelength, emission = plasma.calculate_pcx_chord_emission(loc, cube[0],
        #             w0, mu, Lnu, cube[1], rmax=rmax, nr=101, nlambda=2000,
        #             R_outer=r_anode)
        #     model_signal[idx, :] = models.general_model(rr, L, d, F, wavelength, emission)
        #     model_signal[idx, :] *= cube[image_idx+2]

        # fig, ax = plt.subplots()
        # ax.hist(cubes[:, 1])
        # plt.show()


        fig, ax = plt.subplots()
        ax.plot(rr, ss, color='C0', zorder=10)
        for model_sig in model_signal:
            ax.plot(rr, model_sig, 'C1', zorder=1)
        ax.plot(rr, model_signal.mean(axis=0), 'C2', zorder=11)
        ax.set_xlabel("R (px)")
        ax.set_ylabel("Counts")
        plt.show()
    else:
        pymultinest.run(loglikelihood, prior, n_params,
                importance_nested_sampling=False, resume=resume,
                verbose=True, sampling_efficiency='model',
                n_live_points=150,
                outputfiles_basename=os.path.join(output_folder, 'Ti_chord_')
                )

if __name__ == "__main__":
    # This is for TESTING
    output_folder = None
    calib_posterior = np.random.multivariate_normal([74.6/0.00586, 0.8836181, 21.0], np.diag([0.2/0.00586, 1e-5, 0.2]), size=100)
    print(calib_posterior.shape)
    fname = "/home/milhone/fabry_perot/bin/2019_11_06/scan_2019_11_06_v2.csv"
    #info_df = pd.read_csv(fname, delimiter=',')

    #print(info_dict)
    #print(info_dict['impact_factor'])
    #print(info_df)
    #print(info_df.columns)
    #argon_multi_image_solver_fixed_Lnu(output_folder, calib_posterior, info_dict, image_index=4, resume=True, test_plot=False)
    argon_multi_image_solver_fixed_Lnu(output_folder, calib_posterior, fname, image_index=4, resume=True, test_plot=True)
