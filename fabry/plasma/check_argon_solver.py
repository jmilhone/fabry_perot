from __future__ import print_function, division
import concurrent.futures
import pymultinest
import numpy as np
from ..core import models
from ..tools import plotting, file_io
from . import plasma
import matplotlib.pyplot as plt
import os.path as path


def calculate_profile_model(r, L, d, F, cube, impact_factor, w0, mu, nr, nlambda, Lne, R_outer, r_max):
    w, spec = plasma.calculate_pcx_chord_emission(impact_factor, cube[0], w0, mu,
            cube[3], cube[2], nr=nr, nlambda=nlambda, Lne=Lne, R_outer=R_outer, rmax=r_max)

    vals = cube[1] * models.general_model(r, L, d, F, w, spec)

    return vals



w0 = 487.98634
mu = 39.948


def check_profile_solver(output_folder, Fpost, Lpost, dpost):

    data_filename = path.join(output_folder, "argon_input.h5")
    data = file_io.h5_2_dict(data_filename)

    ix = data['fit_ix']['0'][0:-1:2]
    r = data['r'][ix]
    sig = data['sig'][ix]
    error = data['sig_sd'][ix]

    nL = len(Lpost)
    nF = len(Fpost)

    n_params = 4

    # useful things that I used to make synthetic image but not part of the solver
    impact_factor = 33.0
    nr = 400
    nlambda = 2000
    Lne = 4.0
    R_outer = 35.0
    r_max = 42.0

    # Let's grab the posterior results
    analyzer = pymultinest.Analyzer(n_params, outputfiles_basename=path.join(output_folder, "Ti_profileV_"))
    post = analyzer.get_equal_weighted_posterior()
    nrows, ncols = post.shape

    # Let's sample from all of the posterior results
    n_samples = 3000
    idx_L = np.random.choice(nL, size=n_samples)
    # idx_F = np.random.choice(nF, size=n_samples)

    L_sampled = Lpost[idx_L]
    d_sampled = dpost[idx_L]
    # F_sampled = Fpost[idx_F]
    F_sampled = Fpost[idx_L]

    idx_cube = np.random.choice(nrows, size=n_samples)
    cubes = post[idx_cube, :]

    # test1 = calculate_profile_model(r, 150.0/0.004, 0.88, 20.7, [0.52, 1.0, 4000, 20.0], impact_factor, w0, mu, nr, nlambda, Lne, R_outer, r_max)
    # test2 = calculate_profile_model(r, 150.0/0.004, 0.88, 20.7, [0.5, 1.0, 4000, 20.0], impact_factor, w0, mu, nr, nlambda, Lne, R_outer, r_max)
    # fig, ax = plt.subplots()
    # ax.plot(r, error / sig * 100)
    # plt.show()

    # fig, ax = plt.subplots()
    # ax.plot(r, test1 / test1.max(), 'C0', label='0.52, 4.0')
    # ax.plot(r, test2 / test2.max(), 'C1', label='0.5, 4.0')
    # ax.legend()
    # plt.show()

    vals = np.zeros((n_samples, len(r)))
    # Calculate the model values for all of the samples 
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # submit all of the calcultions to the executor and map to their index
        futures_map = dict()
        for j, (L, d, F, cube) in enumerate(zip(L_sampled, d_sampled, F_sampled, cubes)):
            #new_cube = [0.5, cube[1], 4000, cube[3]]
            fut = executor.submit(calculate_profile_model, r, L, d, F, cube, impact_factor, w0, mu, nr, nlambda, Lne, R_outer, r_max)
            #fut = executor.submit(calculate_profile_model, r, L, d, F, cube, impact_factor, w0, mu, nr, nlambda, Lne, R_outer, r_max)
            futures_map[fut] = j

        # grab the results from completed futures and store in vals
        for future in concurrent.futures.as_completed(futures_map):
            index = futures_map[future]
            vals[index, :] = future.result()
    val_mean = np.nanmean(vals, axis=0)
    # val_std = np.nanstd(vals, axis=0)

    fontsize = 8
    figsize = 3.5
    plot_folder = path.join(output_folder, 'Plots')
    file_io.prep_folder(plot_folder)

    # Make Fit and Data comparison plot
    fig, ax = plt.subplots(figsize=(figsize, figsize/1.618))
    ax.errorbar(r, sig, yerr=error, label='Data', color='C1', errorevery=50, lw=1)
    print(r.shape)
    ax.fill_between(r, np.min(vals, axis=0), np.max(vals, axis=0), color='C0', alpha=0.6, label='Model', hatch='//')
    #ax.plot(r, val_mean, color='C0', label='Model')
    ax.legend(fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax.set_xlabel("R (px)", fontsize=fontsize)
    ax.set_ylabel("Counts", fontsize=fontsize)

    axL = fig.add_axes([0.18+0.1, 0.4+0.1, 0.20, 0.20])
    axR = fig.add_axes([0.63+0.1, 0.4+0.1, 0.20, 0.20])

    axL.tick_params(labelsize=fontsize)
    axR.tick_params(labelsize=fontsize)
    

    i1 = np.abs(r-707.5).argmin()
    i2 = np.abs(r-712.5).argmin()
    i3 = np.abs(r-760).argmin()
    i4 = np.abs(r-765).argmin()

    iL = slice(i1, i2)
    iR = slice(i3, i4)

    axL.errorbar(r[iL], sig[iL], yerr=error[iL], label='Data', color='C1')
    axL.fill_between(r[iL], np.min(vals[:, iL], axis=0), np.max(vals[:, iL], axis=0), color='C0', alpha=0.6, hatch='//')

    axR.errorbar(r[iR], sig[iR], yerr=error[iR], label='Data', color='C1')
    axR.fill_between(r[iR], np.min(vals[:, iR], axis=0), np.max(vals[:, iR], axis=0), color='C0', alpha=0.6, hatch='//')

    #axL.set_ylim(0, 400)
    #axR.set_ylim(0, 400)

    fig.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.18)
    fig.savefig(path.join(plot_folder, 'argon_plasma_fit_vel_profile.pdf'))
    #plt.close(fig)
    plt.show(block=True)

    # Make a plot of the log likelihood
    # chisq = (sig - val_mean)**2 / error**2
    chisq = (sig - vals)**2 / error**2
    #chisq -= np.max(chisq)

    #log_like = -chisq / 2.0
    log_like = chisq / chisq.max()

    fig, ax = plt.subplots()
    # ax.plot(r, log_like)
    ax.fill_between(r, np.min(log_like, axis=0), np.max(log_like, axis=0))
    plt.show()


    # Make Marginal Posterior Plots
    xlabels = ['$T_i$ (eV)', 'A (Counts)', 'V (km/s)', '$L_{\\nu}$ (cm)']
    ylabels = ['$P(T_i) \Delta T_i$', '$P(A) \Delta A$', '$P(V) \Delta V$', '$P(L_{\\nu}) \Delta L_{\\nu}$']
    fnames = ["{}_marginal_vel_profile.png".format(x) for x in ('Ti', 'A', 'V', "Lnu")]
    factors = [1.0, 1.0, 1e-3, 1.0]

    for idx, (xlabel, ylabel, factor) in enumerate(zip(xlabels, ylabels, factors)):
        fig, ax = plt.subplots(figsize=(figsize, figsize/1.618))
        plotting.my_hist(ax, post[:, idx]*factor, bins='auto')
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        fig.savefig(path.join(plot_folder, fnames[idx]))
        # plt.close(fig)
        plt.show(block=True)

    for i, xlabel in enumerate(xlabels):
        for j in range(i+1, 4):
            fig, ax = plt.subplots()
            cb = plotting.my_hist2d(ax, post[:, i]*factors[i], post[:, j]*factors[j])  
            cb.ax.tick_params(labelsize=fontsize)
            ax.set_xlabel(xlabels[j], fontsize=fontsize)
            ax.set_ylabel(xlabel, fontsize=fontsize)
            ax.tick_params(labelsize=fontsize)
            plt.show(block=False)
    plt.show()



