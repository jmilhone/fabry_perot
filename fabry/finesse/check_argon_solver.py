from __future__ import print_function, division, absolute_import
import functools
import concurrent.futures
import pymultinest
import numpy as np
from ..core.models import forward_model, peak_calculator, offset_forward_model
from ..tools import plotting
import json
import argparse
import h5py
import matplotlib.pyplot as plt # For testing purposes only
from ..tools import file_io as io
import random
from os.path import abspath, join

w0 = (487.873302, 487.98634)
mu = (232.03806, 39.948)

def check_solver(finesse_folder, Lpost, dpost):
    """
    Verfies the output of the point spread function solver.  Makes many plots.

    Args:
        finesse_folder (str): path to folder containing output files from the finesse solver
        Lpost (np.ndarray): array of equally weighted marginal posterior L values
        dpost(np.ndarray): array of equally weighted marginal posterior d values
    """
    prior_filename = join(finesse_folder, 'finesse_prior_info.json')
    data_filename = join(finesse_folder, 'finesse_input.h5')

    with open(prior_filename, 'r') as infile:
        prior = json.load(infile, parse_float=np.float64)

    data = io.h5_2_dict(data_filename)

    ix = data['fit_ix']
    r = data['r'][ix]
    print(len(r))
    sig = data['sig'][ix]
    error = data['sig_sd'][ix]

    w_extra = prior['w_extra']

    n_params = 4 + len(w_extra)

    analyzer = pymultinest.Analyzer(n_params, outputfiles_basename=join(finesse_folder, "finesse_"))
    post = analyzer.get_equal_weighted_posterior()

    F_post = post[:, 0]
    A_post = post[:, 1]
    Arel_post = post[:, 2]
    Ti_post = post[:, 3]

    # subsample Lpost and dpost
    nL = 50 
    Ld_choice = np.random.choice(len(Lpost), size=nL)
    L_vals = Lpost[Ld_choice]
    d_vals = dpost[Ld_choice]

    # subsample from finesse posterior
    nF = 150
    post_choice = np.random.choice(len(F_post), size=nF)
    F_vals = post[post_choice, 0]
    A_vals = post[post_choice, 1]
    Arel_vals = post[post_choice, 2]
    Ti_vals = post[post_choice, 3]

    output = np.zeros((len(r), nL, nF))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        idx = 0
        output_index_map = {}
        futures_map = {}
        for i, (L, d) in enumerate(zip(L_vals, d_vals)):
            for j, (F, A, Arel, Ti) in enumerate(zip(F_vals, A_vals, Arel_vals, Ti_vals)):
                fut = executor.submit(offset_forward_model, r, L, d, F, w0, mu, [A*Arel, A], 
                        [0.025, Ti], [0.0, 0.0])
                #fut = executor.submit(int, idx)
                futures_map[fut] = idx
                output_index_map[idx] = (i, j)
                idx += 1

        # for k, v in futures_map.items():
        #     print(k ,v)

        for future in concurrent.futures.as_completed(futures_map):
            idx = futures_map[future]
            i, j = output_index_map[idx]

            output[:, i, j] = future.result()

    output.shape = (len(r), nL*nF)
    output_mean = np.mean(output, axis=1)
    output_max = np.max(output, axis=1)
    output_min = np.min(output, axis=1)

    #for ii in range(nL*nF):
    #    ax.plot(r, output[:, ii], 'C0')
    fig, ax = plt.subplots()
    ax.fill_between(r, output_min, output_max, color='C0', alpha=0.5)
    ax.plot(r, output_mean, 'C0')
    ax.plot(r, sig, 'C1')
    plt.show()

    fig, ax = plt.subplots()
    ax.hist(peak_calculator(L_vals, d_vals, w0[0], 0), bins='auto')
    plt.show()
    # print(L, d, F, A, Arel, Ti)
    # vals = forward_model(r, L, d, mode[0], w, mass, amps, Ti, V, sm_ang=False, nlambda=2000)
    # mode_stats = analyzer.get_mode_stats()

    # # Find the mode with the lowest local log-evidence
    # minval = np.inf
    # ix = 0
    # for idx, stat in enumerate(mode_stats['modes']):
    #     if stat['local log-evidence'] < minval:
    #         print('im here')
    #         ix = idx
    #         minval = stat['local log-evidence']

    # mode = mode_stats['modes'][ix]['mean']

    # amps = np.zeros(2 + len(w_extra))
    # amps[0] = mode[2]
    # amps[1] = 1.0
    # amps[2:] = mode[4:]
    # amps *= mode[1]
    # amps = list(amps)

    # w = [x for x in w0]
    # w += w_extra

    # mass = [x for x in mu]
    # mass += [mu[0] for _ in w_extra]

    # V = [0.0 for _ in mass]

    # Ti = [0.025*1000.0/300.0, mode[3]]
    # Ti += [0.025*1000.0/300.0 for _ in w_extra]

    # #L = 0.380173301412519577E+05
    # #d = 0.883628502371783142E+00
    # ichoice = random.choice(range(len(Lpost)))
    # L = Lpost[ichoice]
    # d = dpost[ichoice]
    # vals = forward_model(r, L, d, mode[0], w, mass, amps, Ti, V, sm_ang=False, nlambda=2000)

    # # trying to qmodel offset here
    # vals += mode[1] * 0.15 / (1.0 + mode[0])

    # fig, ax = plt.subplots()
    # ax.errorbar(r, sig, yerr=error, ecolor='C2', color='C0', label='Data')
    # ax.plot(r, vals, color='C1', label='Fit')
    # styles = ['-', '--', '-.']
    # for idx, wavelength in enumerate(w):
    #     n, m = divmod(idx, 10)
    #     if wavelength in [488.185311, 488.32348]:
    #         order = 1
    #     else:
    #         order = 0
    #     ax.axvline(peak_calculator(L, d, wavelength, order), label=str(wavelength), color='C{0:d}'.format(m), linestyle=styles[n])
    # ax.legend()
    # plt.show(block=False)

    # fig, ax = plt.subplots()
    # ax.plot(r, (vals - sig)**2 / error**2)
    # plt.show()

def check_full_solver(finesse_folder):
    """
    Verfies the output of the full point spread function solver.  Makes many plots.

    Args:
        finesse_folder (str): path to folder containing output files from the finesse solver
    """

    data_filename = join(finesse_folder, 'finesse_input.h5')

    data = io.h5_2_dict(data_filename)

    ix = data['fit_ix']
    r = data['r'][ix]
    sig = data['sig'][ix]
    error = data['sig_sd'][ix]


    n_params = 6

    analyzer = pymultinest.Analyzer(n_params, outputfiles_basename=join(finesse_folder, "full_"))
    post = analyzer.get_equal_weighted_posterior()


    nrows, ncols = post.shape 
    print(nrows, ncols)
    print(len(r))
    L = post[:, 0]
    d = post[:, 1]
    F = post[:, 2]
    A = post[:, 3]
    Arel = post[:, 4]
    Ti = post[:, 5]
    amp = [Arel*A, A]
    temps = [0.025, Ti]
    V = [0.0, 0.0]

    npts=nrows
    i = np.random.choice(nrows, size=npts)
    vals = np.zeros((npts, len(r)))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures_map = {}
        for j, idx in enumerate(i):
            fut = executor.submit(offset_forward_model, r, L[idx], d[idx], F[idx], w0, mu,
                    [Arel[idx]*A[idx], A[idx]], [0.025, Ti[idx]], V, nlambda=2000)
            futures_map[fut] = j
        for future in concurrent.futures.as_completed(futures_map):
            idx = futures_map[future]
            vals[idx, :] = future.result()

    val_mean = np.nanmean(vals, axis=0)
    val_std = np.nanstd(vals, axis=0)
    #vals = offset_forward_model(r, L, d, F, list(w0), list(mu), amp, temps, V, nlambda=2000)

    fig, ax = plt.subplots()
    ax.errorbar(data['r'], data['sig'], yerr=data['sig_sd'], label='data', color='C1')
    #ax.fill_between(r, val_mean-val_std, val_mean+val_std, color='C0', alpha=0.6)
    ax.fill_between(r, np.min(vals, axis=0), np.max(vals, axis=0), color='C0', alpha=0.6)
    ax.plot(r, val_mean, 'C0', label='fit')
    ax.legend()
    plt.show()

    pk = peak_calculator(L, d, w0[1], order=0)

    labels = ["L", "d", "F", "A", "Arel", "Ti"]
    fac = [0.004, 1.0, 1.0, 1.0, 1.0, 1.0]
    bins = ['auto', 50, 'auto', 'auto', 'auto', 'auto']
    for idx, label in enumerate(labels):
        print(label, ": ", np.nanmean(post[:, idx]*fac[idx]))
        fig, ax = plt.subplots()
        plotting.my_hist(ax, post[:, idx]*fac[idx], bins=bins[idx])
        ax.set_xlabel(label)
        plt.show(block=False)
    plt.show(block=True)
