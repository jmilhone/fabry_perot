from __future__ import print_function, division

import os.path

import numpy as np
import pymultinest

import fabry.core.models as models
import fabry.distributions.distribution as dist
# Go back and look into relative imports!
# from fabry.tools import file_io
import fabry.plasma.plasma as plasma

w0 = 487.98634
mu = 39.948

# For testing purposes, dont' have a calibration to use
L = 74.6 / 0.00586
d = 0.8836181
F = 20.9


def pcx_vfd_model(r, b, Ti_params, V_params, ne_params, delta_d, amp, L, d, F, rmax=40.0):
    """

    :param r:
    :param b:
    :param Ti_params:
    :param V_params:
    :param ne_params:
    :param delta_d:
    :param amp:
    :param L:
    :param d:
    :param F:
    :param rmax:
    :return:
    """
    w_model, radiance_model = plasma.vfd_chord(b, Ti_params, V_params, ne_params,
                                               rmax=rmax, w0=w0, mu=mu, nr=101,
                                               nlambda=2048, test_plot=False)

    model_signal = amp * models.general_model(r, L, d + delta_d, F, w_model, radiance_model)
    return model_signal


def calculate_delta_d(velocity_offset):
    """Calculates delta d (in mm) for a velocity offset (m/s)

    :param float velocity_offset: offset velocity in n/s
    :return: delta d in mm
    :rtype: float
    """
    return 2.94734982e-09 * velocity_offset


def argon_multi_image_solver_fixed_Lnu(output_folder, calib_posterior,
                                       image_data, resume=True, test_plot=False):
    """Runs MultiNest solver for multiple chords in a PCX Ar plasma with a fixed Lnu.

    :param str output_folder: output folder to store multinest files
    :param np.ndarray calib_posterior: calibration equally weighted posterior samples
    :param tuple image_data: Tuple containing (r, sig, std, impact_factor) for the images to be analyzed
    :param bool resume: Resume multinest solver if True, default=True
    :param bool test_plot: Runs test plot sampling from prior instead of running solver, default=False
    """

    def prior(cube, ndim, nparams):
        """Transforms prior hypercube to physical units
        """
        # cube[0] = dist.UniformPrior(cube[0], 0.025, 3.0)
        # cube[1] = dist.UniformPrior(cube[1], -5000.0, 5000.0)
        cube[0] = dist.UniformPrior(cube[0], 0.025, 3.0)
        cube[1] = dist.UniformPrior(cube[1], 0.025, 3.0)
        cube[2] = dist.UniformPrior(cube[2], -5000.0, 5000.0)

        # I want to use a less broad prior for the amplitudes, but not allow negative amplitude
        # for i in range(2, nparams):
        for i in range(3, nparams):
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


        # For each image, calculate the model signal and update the total chisquared
        for i, (loc, rr, ss, std, v_offset) in enumerate(zip(impact_factors, r, sig, sd, vel_offsets)):
            # prepare parameters for the model function
            Ti_params = [cube[0], cube[1]]
            V_params = [cube[2], Lnu, r_anode]
            ne_params = [r_anode, ]
            delta_d = calculate_delta_d(v_offset)
            amp = cube[i + 3]

            mod_signal = pcx_vfd_model(rr, loc, Ti_params, V_params, ne_params,
                                       delta_d, amp, L, d, F, rmax=rmax)
            chisq += np.sum((mod_signal - ss) ** 2 / std ** 2)

        return -chisq / 2.0

    def create_test_plot(n_samples=500, image_idx=0):
        import matplotlib.pyplot as plt
        # sample from prior and see if things make sense
        # L = 74.6 / px_size
        # d = 0.8836181
        # F = 21.0
        cubes = np.random.rand(n_samples, n_params)

        # grab image data to compare prior samples to
        image_idx = 2
        loc = impact_factors[image_idx]
        rr = r[image_idx]
        ss = sig[image_idx]
        std = sd[image_idx]

        model_signal = np.zeros((n_samples, len(rr)))

        for idx, cube in enumerate(cubes):
            # transform unit hypercube to physical units
            prior(cube, None, n_params)

            Ti_params = [cube[0], cube[1]]
            V_params = [cube[2], rmax, r_anode]
            ne_params = [r_anode, ]
            delta_d = calculate_delta_d(vel_offsets[image_idx])
            amp = cube[image_idx + 3]
            model_signal[idx, :] = pcx_vfd_model(rr, loc, Ti_params, V_params,
                                                 ne_params, delta_d, amp, L, d, F,
                                                 rmax=rmax)

        # fig, ax = plt.subplots()
        # ax.hist(cubes[:, 1])
        # plt.show()

        fig, ax = plt.subplots()
        ax.fill_between(rr, ss - std, ss + std, color='C0', alpha=0.5)
        ax.plot(rr, ss, color='C0', zorder=10)
        for model_sig in model_signal:
            ax.plot(rr, model_sig, 'C1', zorder=1)
        ax.plot(rr, model_signal.mean(axis=0), 'C2', zorder=11)
        ax.set_xlabel("R (px)")
        ax.set_ylabel("Counts")
        plt.show()

    #### Constants
    px_size = 0.00586
    Lnu = 50.0
    rmax = 40.0
    r_anode = 32.0
    # n_calib_points = calib_posterior.shape[0]

    # unpact image_data
    r = image_data[0]
    sig = image_data[1]
    sd = image_data[2]
    impact_factors = image_data[3]
    n_images = len(r)

    print("Code in the velocity offsets!")
    vel_offsets = [0.0 for _ in range(n_images)]

    # params are Ti, Vmax, +amplitudes for each image
    # n_params = 2 + n_images

    # params are Ti_inner, Ti_outer, Vmax, +amplitudes for each image
    n_params = 3 + n_images

    if test_plot:
        create_test_plot(n_samples=500, image_idx=0)
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
    cal_posterior = np.random.multivariate_normal([74.6 / 0.00586, 0.8836181, 21.0],
                                                    np.diag([0.2 / 0.00586, 1e-5, 0.2]), size=100)
    print(cal_posterior.shape)
    fname = "/home/milhone/fabry_perot/bin/2019_11_06/scan_2019_11_06_v2.csv"
    # info_df = pd.read_csv(fname, delimiter=',')

    # print(info_dict)
    # print(info_dict['impact_factor'])
    # print(info_df)
    # print(info_df.columns)
    # argon_multi_image_solver_fixed_Lnu(output_folder, cal_posterior, info_dict, image_index=4, resume=True, test_plot=False)
    argon_multi_image_solver_fixed_Lnu(output_folder, cal_posterior, fname, resume=True,
                                       test_plot=True)
