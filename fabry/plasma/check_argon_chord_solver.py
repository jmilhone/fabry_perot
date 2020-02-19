from __future__ import print_function, division

import concurrent.futures
import os.path as path

import matplotlib.pyplot as plt
import numpy as np
import pymultinest
from matplotlib import rcParams

import fabry.plasma.argon_chord_solver as acs
from fabry.tools import plotting, file_io

rcParams['font.size'] = 16
rcParams['figure.figsize'] = (8, 6)

w0 = 487.98634
mu = 39.94


def check_const_Lnu_solver(output_folder, calib_posterior, image_data, n_samples=200):
    plot_folder = path.join(output_folder, 'Plots')
    file_io.prep_folder(plot_folder)

    # Hard coded calibration!
    hard_coded_calib = True
    if hard_coded_calib:
        print('Using a hard coded calibration')
        from fabry.plasma.argon_chord_solver import L, d, F
        print("L, d, F = ", L, d, F)
    else:
        idx_calib = np.random.choice(calib_posterior.shape[0], size=n_samples)
        L = calib_posterior[idx_calib, 0]
        d = calib_posterior[idx_calib, 1]
        F = calib_posterior[idx_calib, 2]

    px_size = 0.00586
    Lnu = 100.0
    rmax = 40.0
    r_anode = 32.0

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

    # params are Ti inner, Ti_outer, Vmax, +amplitudes for each image
    n_params = 3 + n_images

    analyzer = pymultinest.Analyzer(n_params, outputfiles_basename=path.join(output_folder, "Ti_chord_"))
    post = analyzer.get_equal_weighted_posterior()
    # print("overwriting velocity in post!")
    # post[:, 1] -= -200
    post_means = np.mean(post, axis=0)
    print(post_means)
    nrows, ncols = post.shape

    idx_cube = np.random.choice(nrows, size=n_samples)
    cubes = post[idx_cube, :]


    # allocate memory for posterior values
    values = []
    for rr in r:
        values.append(np.zeros((n_samples, len(rr))))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures_map = dict()
        # for j, (L, d, F, cube) in enumerate(zip(L_sampled, d_sampled, F_sampled, cubes)):
        for j, cube in enumerate(cubes):
            fut = executor.submit(calculate_multi_models, r, impact_factors, L, d, F, Lnu, cube)
            futures_map[fut] = j

        for future in concurrent.futures.as_completed(futures_map):
            index = futures_map[future]
            output = future.result()
            for idx, out in enumerate(output):
                values[idx][index, :] = out

    #################################################### 
    ## Calculate percentiles for the posterior values ##
    #################################################### 
    levels = [68, 95, 99]
    val_percentiles = list()
    val_mean = list()
    for value in values:
        val_percentiles.append(calculate_percentile_ranges(value, levels))
        val_mean.append(np.mean(value, axis=0))

    ############## 
    ## Fit Plot ##
    ############## 
    alphas = [0.5, 0.3, 0.1]
    for i, (rr, ss, std, val_percentile) in enumerate(zip(r, sig, sd, val_percentiles)):
        fig, ax = plt.subplots()
        for alpha, level in zip(alphas, levels):
            ax.fill_between(rr, val_percentile[level][0], val_percentile[level][1], color='k', alpha=alpha)
        #ax.fill_between(rr, values[i].min(axis=0), values[i].max(axis=0), color='k', alpha=0.3)
        ax.plot(rr, val_mean[i], color='k', label='Fit')
        ax.errorbar(rr, ss, yerr=std, color='C1', errorevery=5, label='Data')
        impact_factor = impact_factors[i]

        title_string = f"R = {impact_factor:2.1f} cm"
        ax.set_title(title_string)
        ax.legend()
        plt.show()

    ######################## 
    ## Histogram Plotting ##
    ######################## 
    xlabels = ['$T_i$ (eV)', '$V_{pk}$ (km/s)']
    xlabels = ['$T_i$ (eV)', '$V_{pk}$ (km/s)']
    xlabels = ['$T_i$ (eV)', '$V_{pk}$ (km/s)']
    xlabels = ['$T_i$ (eV)', '$V_{pk}$ (km/s)']
    xlabels = ['$T_i$ (eV)', '$V_{pk}$ (km/s)']
    ylabels = ['$P(T_i) \Delta T_i$ (%)', '$P(V_{pk}) \Delta V_{pk}$ (%)']
    factors = [1.0, 1e-3]

    # add in the amplitudes
    xlabels += ["$A_{0}$".format(i) for i in range(n_images)]
    ylabels += ["$P(A_{0}) \Delta A_{0}$ (%)".format(i) for i in range(n_images)]
    factors += [1.0 for _ in range(n_images)]

    make_histograms(post, xlabels, ylabels, factors, None, save=False)

    ##################################################### 
    ## 2D Marginal Posterior Distribution for Ti and V ##
    ##################################################### 
    fig, ax = plt.subplots()
    cb = plotting.my_hist2d(ax, post[:, 0]*factors[0], post[:, 1]*factors[1])
    #cb.ax.tick_params(labelsize=16)
    ax.set_xlabel(xlabels[1])
    ax.set_ylabel(xlabels[0])
    plt.show()


def make_histograms(posterior, xlabels, ylabels, factors, fnames, save=False, plot_folder=None):

    for idx, (xlabel, ylabel, factor) in enumerate(zip(xlabels, ylabels, factors)):
        fig, ax = plt.subplots()
        plotting.my_hist(ax, posterior[:, idx]*factor, bins='auto')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        if save and plot_folder is not None:
            filepath = path.join(plot_folder, fnames[idx])
            fig.savefig(filepath, transparent=True)

        plt.show()

def calculate_multi_models(r, impact_factors, vel_offsets, L, d, F, Lnu, cube, nlambda=2000,
                           nr=101, rmax=40.0, r_anode=32.0):
    values = list()
    for idx, (loc, v_offset) in enumerate(zip(impact_factors, vel_offsets)):
        # wavelength, emission = plasma.calculate_pcx_chord_emission(loc, cube[0],
        #         w0, mu, Lnu, cube[1], rmax=rmax, nr=101, nlambda=2000,
        #         R_outer=r_anode)
        # model_signal = models.general_model(r[idx], L, d, F, wavelength, emission)
        # model_signal = cube[idx+2] * model_signal

        # prepare parameters for the model function
        Ti_params = [cube[0], cube[1]]
        V_params = [cube[2], Lnu, r_anode]
        ne_params = [r_anode, ]
        delta_d = acs.calculate_delta_d(v_offset)
        amp = cube[idx + 3]

        model_signal = acs.pcx_vfd_model(r[idx], loc, Ti_params, V_params, ne_params,
                                                                     delta_d, amp, L, d, F, rmax=rmax)
        values.append(model_signal)

    return values

def calculate_percentile_ranges(data, levels):
    lower_levels = [50.0 - x/2.0 for x in levels]
    upper_levels = [50.0 + x/2.0 for x in levels]

    percentiles = {}
    for level, LL, UL in zip(levels, lower_levels, upper_levels):
        lower = np.percentile(data, LL, axis=0)
        upper = np.percentile(data, UL, axis=0)
        percentiles[level] = (lower, upper)

    return percentiles








