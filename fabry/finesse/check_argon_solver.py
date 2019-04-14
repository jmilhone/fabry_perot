from __future__ import print_function, division, absolute_import
import functools
import concurrent.futures
import pymultinest
import numpy as np
from ..core.models import forward_model, peak_calculator, offset_forward_model
from ..tools import plotting
import json
import h5py
import matplotlib.pyplot as plt # For testing purposes only
from ..tools import file_io as io
import random
from os.path import abspath, join

w0 = (487.873302, 487.98634, 487.800942)
mu = (232.03806, 39.948, 232.03806)

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

    ix = data['fit_ix']['0']
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
                fut = executor.submit(forward_model, r, L, d, F, w0, mu, [A*Arel, A], 
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
    test = forward_model(r, 150.0/0.004, 0.88, 20.7, w0, mu, [0.5*1000.0, 1000.0], [0.025, 0.1], [0.0, 0.0])
    fig, ax = plt.subplots()
    ax.fill_between(r, output_min, output_max, color='C0', alpha=0.5)
    ax.plot(r, output_mean, 'C0')
    ax.plot(r, sig, 'C1')
    #ax.plot(r, test, 'C2', labe
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

    ix = data['fit_ix']['0']
    r = data['r'][ix]
    sig = data['sig'][ix]
    error = data['sig_sd'][ix]


    #n_params = 8
    n_params = 7
    #n_params = 6

    analyzer = pymultinest.Analyzer(n_params, outputfiles_basename=join(finesse_folder, "full_"))
    post = analyzer.get_equal_weighted_posterior()

    nrows, ncols = post.shape
    print('Posterior Shape: ({0},{1})'.format(nrows, ncols))

    # Unpack posterior into relevant variable names
    L = post[:, 0]#*3
    d = post[:, 1]
    F = post[:, 2]
    A = post[:, 3]
    Arel = post[:, 4]
    Ti = post[:, 5]
    #print('****************************')
    #print('offset set to zero for Cal6')
    #print('****************************')
    #offset = post[:, 6]
    #offset = np.zeros_like(Ti)
    Brel = post[:, 6]
    #Brel = post[:, 7]
    amp = [Arel*A, A]
    temps = [0.025, Ti]
    V = [0.0, 0.0, 0.0]

    npts=nrows
    i = np.random.choice(nrows, size=npts)
    vals = np.zeros((npts, len(r)))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures_map = {}
        for j, idx in enumerate(i):
            fut = executor.submit(forward_model_wrapper, r, L[idx], d[idx],
                    F[idx], w0, mu,[Arel[idx]*A[idx], A[idx], Brel[idx]*A[idx]],
                    [0.085, Ti[idx], 0.085], V, 0.0, nlambda=2000)
            #fut = executor.submit(forward_model_wrapper, r, L[idx], d[idx],
            #        F[idx], w0, mu,[Arel[idx]*A[idx], A[idx]],
            #        [0.085, Ti[idx]], V, 0.0, nlambda=2000)#offset[idx], nlambda=2000)
            futures_map[fut] = j
        for future in concurrent.futures.as_completed(futures_map):
            idx = futures_map[future]
            vals[idx, :] = future.result()

    val_mean = np.nanmean(vals, axis=0)
    val_std = np.nanstd(vals, axis=0)

    output_fit = {'r': r, 'sig': val_mean}
    io.dict_2_h5(join(finesse_folder, "fit_data.h5"), output_fit)

    plot_folder = join(finesse_folder, 'Plots')
    io.prep_folder(plot_folder)

    # idx = np.random.choice(post.shape[0])
    # L = post[idx, 0]#*3
    # d = post[idx, 1]
    # F = post[idx, 2]
    # A = post[idx, 3]
    # Arel = post[idx, 4]
    # Ti = post[idx, 5]
    # offset = post[idx, 6]
    # #offset = 0.0
    # d += 487.98634e-6 / 2.0 # Move d by one order
    
    #test_model = forward_model_wrapper(r, L, d, F, w0, mu, [Arel*A, A], [0.085, Ti], V, offset, nlambda=2000)

    ###################### 
    ###  Plot Results  ###
    ###################### 

    # Fit Plot
    fig, ax = plt.subplots(figsize=(12,8))
    ax.errorbar(data['r'], data['sig'], yerr=data['sig_sd'], label='data', color='C1')
    ax.fill_between(r, np.min(vals, axis=0), np.max(vals, axis=0), color='C0', alpha=0.6)
    ax.plot(r, val_mean, 'C0', label='fit')
   # ax.plot(r, test_model, 'C4', label='d+1order')
    ax.legend(fontsize=16)
    ax.set_xlabel("R (px)", fontsize=16)
    ax.set_ylabel("Counts", fontsize=16)
    ax.tick_params(labelsize=16)
    ax.set_xlim(r.min()*0.9, r.max()*1.1)
    plt.show()
    
    idx = np.abs(r - 640).argmin()
    fig, ax = plt.subplots()
    ax.hist(vals[:, idx], bins='auto')
    plt.show()

    labels = ["L", "d", "F", "A", "Arel", "Ti"]#, 'Offset', "Brel"]
    fac = [0.004, 1.0, 1.0, 1.0, 1.0, 1.0, ]#1.0, 1.0]

    #  Use lambda / 10 for a bin size 
    dbins = (np.max(d)-np.min(d)) / (488e-6/10)
    print(np.min(d), np.max(d))
    dbins = int(dbins)
    print('dbins', dbins)
    #dbins = 100
    if dbins < 2:
        dbins = 50
    print(dbins)
    # 1D Marginal Distributions
    bins = ['auto', dbins, 'auto', 'auto', 'auto', 'auto', ]#'auto', 'auto']
    fnames = ["{}_marginal.png".format(x) for x in ('L', 'd', 'F', 'A', 'Arel', "Ti")]#, 'Offset', 'Brel')]
    for idx, label in enumerate(labels):
        print(label, ": ", np.nanmean(post[:, idx]*fac[idx]))
        fig, ax = plt.subplots()
        plotting.my_hist(ax, post[:, idx]*fac[idx], bins=bins[idx])
        ax.set_xlabel(label)
        fig.savefig(join(plot_folder, fnames[idx]))
        plt.show(block=False)
    plt.show(block=True)

    fontsize=16
    for i, xlabel in enumerate(labels):
        for j in range(i+1, 6):
            fname = "{0}_{1}_joint_marginal_dist.pdf".format(xlabel, labels[j])
            fig, ax = plt.subplots()
            cb = plotting.my_hist2d(ax, post[:, i]*fac[i], post[:, j]*fac[j])  
            cb.ax.tick_params(labelsize=fontsize)
            ax.set_xlabel(labels[j], fontsize=fontsize)
            ax.set_ylabel(xlabel, fontsize=fontsize)
            ax.tick_params(labelsize=fontsize)
            fig.tight_layout()
            fig.savefig(join(plot_folder, fname), transparent=True)
            plt.close(fig)
            #plt.show(block=False)
    #plt.show()


def forward_model_wrapper(r, L, d, F, w0, mu, amp, temps, V, offset, nlambda=2000):
    vals = forward_model(r, L, d, F, w0, mu, amp, temps, V, nlambda=nlambda)
    vals += offset
    return vals
