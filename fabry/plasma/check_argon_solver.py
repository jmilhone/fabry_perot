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

    Lnu = 100.0 * plasma.Lnu(cube[3], cube[0], mu=40, noise=False)
    w, spec = plasma.calculate_pcx_chord_emission(impact_factor, cube[0], w0, mu,
            Lnu, cube[2], nr=nr, nlambda=nlambda, Lne=Lne, R_outer=R_outer, rmax=r_max)

    vals = cube[1] * models.general_model(r, L, d, F, w, spec)

    return vals

def const_model(r, L, d, F, w0, mu, cube, nlambda):
    #vals = models.forward_model(r, L, d, cube[3], w0, mu, cube[1], cube[0], cube[2], nlambda=nlambda)
    vals = models.forward_model(r, L, d, F, w0, mu, cube[1], cube[0], cube[2], nlambda=nlambda)
    #vals += 2.7#cube[3]
    return vals



w0 = 487.98634
mu = 39.948


def check_profile_solver(output_folder, Fpost, Lpost, dpost):
    plot_folder = path.join(output_folder, 'Plots')
    file_io.prep_folder(plot_folder)

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
    impact_factor = 25.0
    nr = 400
    nlambda = 2000
    Lne = 4.0
    R_outer = 35.0
    r_max = 42.0

    # Let's grab the posterior results
    analyzer = pymultinest.Analyzer(n_params, outputfiles_basename=path.join(output_folder, "Ti_profileV_"))
    post = analyzer.get_equal_weighted_posterior()
    nrows, ncols = post.shape

    Lnu = np.zeros(nrows)
    rarr = np.linspace(0.0, 50.0, 1001.0)
    v_profile = np.zeros((nrows, len(rarr)))

    for idx, (nen0, Ti, Vouter) in enumerate(zip(post[:, 3], post[:, 0], post[:, 2])):
        Lnu[idx] = 100.0 * plasma.Lnu(nen0, Ti, mu=40, noise=False)
        v_profile[idx, :] = plasma.pcx_velocity_profile(rarr, Lnu[idx], R_outer, Vouter)

    fig, ax = plt.subplots()
    ax.hist(Lnu)
    plt.show()

    v_perc = calculate_percentile_ranges(v_profile)
    levels = [68, 95, 99]
    alphas = [0.8, 0.5, 0.2]

    figsize = (3.5, 3.5/1.618)
    fontsize = 9
    fig, ax = plt.subplots(figsize=figsize)
    for level, alpha in zip(levels, alphas):
        ax.fill_between(rarr, v_perc[level][0]/1000.0, v_perc[level][1]/1000.0, alpha=alpha, color='C0', label='{0:d}%'.format(level))
    mean_v = np.mean(v_profile, axis=0)
    #ax.fill_between(rarr, np.min(v_profile, axis=0), np.max(v_profile, axis=0), color='C0', alpha=0.5)
    ax.plot(rarr, mean_v/1000.0, 'C1', label='Mean')
    ax.legend()
    ax.tick_params(labelsize=fontsize)
    ax.set_xlabel("R (cm)", fontsize=fontsize)
    ax.set_ylabel("V (km/s)", fontsize=fontsize)
    fig.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.18)
    fig.savefig(path.join(plot_folder, "vel_profile_prediction.pdf"), transparent=True)
    plt.show()
 
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
            #print(impact_factor)
            fut = executor.submit(calculate_profile_model, r, L, d, F, cube, impact_factor, w0, mu, nr, nlambda, Lne, R_outer, r_max)
            #fut = executor.submit(calculate_profile_model, r, L, d, F, cube, impact_factor, w0, mu, nr, nlambda, Lne, R_outer, r_max)
            futures_map[fut] = j

        # grab the results from completed futures and store in vals
        for future in concurrent.futures.as_completed(futures_map):
            index = futures_map[future]
            vals[index, :] = future.result()
    val_mean = np.nanmean(vals, axis=0)
    # val_std = np.nanstd(vals, axis=0)

    imax = np.argmax(val_mean)
    print(r[imax])

    fontsize = 8
    figsize = 3.5

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
    fig.savefig(path.join(plot_folder, 'argon_plasma_fit_vel_profile.pdf'), transparent=True)
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
    xlabels = ['$T_i$ (eV)', 'A (Counts)', '$V_{pk}$ (km/s)', '$L_{\\nu}$ (cm)']
    ylabels = ['$P(T_i) \Delta T_i$', '$P(A) \Delta A$', '$P(V_{pk}) \Delta V_{pk}$', '$P(L_{\\nu}) \Delta L_{\\nu}$']
    # fnames = ["{}_marginal_vel_profile.png".format(x) for x in ('Ti', 'A', 'V', "Lnu")]
    fnames = ["{}_marginal_vel_profile.pdf".format(x) for x in ('Ti', 'A', 'V', "Lnu")]
    factors = [1.0, 1.0, 1e-3, 1.0]
    answers = [0.5, None, 1.0, None,]
    for idx, (xlabel, ylabel, factor) in enumerate(zip(xlabels, ylabels, factors)):
        fig, ax = plt.subplots(figsize=(figsize, figsize/1.618))
        plotting.my_hist(ax, post[:, idx]*factor, bins='auto')
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        if answers[idx] is not None:
            ax.plot(answers[idx], 8.0, '*', ms=10, color='C3')
        fig.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.18)
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
            fig.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.18)
            plt.show(block=False)
    plt.show()



def check_constV_solver(output_folder, Fpost, Lpost, dpost, use_plasma_F=False):
    data_filename = path.join(output_folder, "argon_input.h5")
    data = file_io.h5_2_dict(data_filename)

    ix = data['fit_ix']['0']#[0:-1:2]
    r = data['r']#[ix]
    sig = data['sig']#[ix]
    error = data['sig_sd']#[ix]

    nL = len(Lpost)
    nF = len(Fpost)
    print('F mean', np.mean(Fpost))
    if use_plasma_F:
        n_params = 4
    else:
        n_params = 4
    #n_params = 3
    n_params = 4
    # useful things that I used to make synthetic image but not part of the solver
    nlambda = 2000

    # Let's grab the posterior results
    analyzer = pymultinest.Analyzer(n_params, outputfiles_basename=path.join(output_folder, "Ti_constV_"))
    post = analyzer.get_equal_weighted_posterior()
    nrows, ncols = post.shape

    # Let's sample from all of the posterior results
    n_samples = 1000
    idx_L = np.random.choice(nL, size=n_samples)
    # idx_F = np.random.choice(nF, size=n_samples)

    L_sampled = Lpost[idx_L]#*3
    print(L_sampled.min(), L_sampled.max())
    d_sampled = dpost[idx_L]
    F_sampled = Fpost[idx_L]
    #F_sampled = Fpost[idx_F]

    idx_cube = np.random.choice(nrows, size=n_samples)
    cubes = post[idx_cube, :]

    if False:#use_plasma_F:
        #F_folder = "/home/milhone/Research/python_FabryPerot/Data/PCX/2158/"
        #Fpost = np.loadtxt(path.join(F_folder, 'Ti_constV_post_equal_weights.dat'), ndmin=2)[:,3] 
        #F_folder = "/home/milhone/Research/python_FabryPerot/Data/2018_10_28/"
        #Fpost = np.loadtxt(path.join(F_folder, 'full_post_equal_weights.dat'), ndmin=2)[:,2]
        nF = len(Fpost)
        idx_F = np.random.choice(nF, size=n_samples)
        F_sampled = Fpost[idx_F]
    else:
        pass
        #F_sampled = cubes[:, 3]


    rr = np.linspace(400, 850, 1000)
    vals = np.zeros((n_samples, len(rr)))
    # Calculate the model values for all of the samples 
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # submit all of the calcultions to the executor and map to their index
        futures_map = dict()
        for j, (L, d, F, cube) in enumerate(zip(L_sampled, d_sampled, F_sampled, cubes)):
            # L = 0.124073186132456485E+05
            # d = 0.876551993778090566E+00
            # F = 0.205637295793309107E+02
        #for j, (L, d, cube) in enumerate(zip(L_sampled, d_sampled, cubes)):
            fut = executor.submit(const_model, rr, L, d, 22.7, w0, mu, cube, nlambda)
            #fut = executor.submit(models.forward_model, r, L, d, F, w0, mu, cube[1], cube[0], cube[2], nlambda=nlambda)
            futures_map[fut] = j

        # grab the results from completed futures and store in vals
        for future in concurrent.futures.as_completed(futures_map):
            index = futures_map[future]
            vals[index, :] = future.result() + cubes[index, 3]
    val_mean = np.nanmean(vals, axis=0)
    imax = np.argmax(val_mean)
    #print(r[imax])

    # val_std = np.nanstd(vals, axis=0)
    mytest = const_model(r, 0.125379300991324253E+05, 0.879962652610075002E+00,
            0.0, w0, mu, [0.0001, 0.43e3, -14000.0, 27.0], nlambda)
    fontsize = 8
    figsize = 3.5
    plot_folder = path.join(output_folder, 'Plots')
    file_io.prep_folder(plot_folder)

    # Make Fit and Data comparison plot
    fig, ax = plt.subplots(figsize=(figsize, figsize/1.618))
    ax.errorbar(r, sig, yerr=error, label='Data', color='C1', errorevery=10, lw=1)
    print(r.shape)
    ax.fill_between(rr, np.min(vals, axis=0), np.max(vals, axis=0), color='C0', alpha=0.6, label='Model', hatch='//')
    #ax.plot(r, val_mean, color='C0', label='Model')
    #ax.plot(r, mytest, 'k')
    ax.legend(fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax.set_xlabel("R (px)", fontsize=fontsize)
    ax.set_ylabel("Counts", fontsize=fontsize)

    fig.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.18)
    fig.savefig(path.join(plot_folder, 'argon_plasma_constV.pdf'), transparent=True)
    #plt.close(fig)
    plt.show(block=True)

    # Make Marginal Posterior Plots
    xlabels = ['$T_i$ (eV)', 'A (Counts)', 'V (km/s)']
    ylabels = ['$P(T_i) \Delta T_i$', '$P(A) \Delta A$', '$P(V) \Delta V$']
    fnames = ["{}_marginal_vel_profile.png".format(x) for x in ('Ti', 'A', 'V')]
    #fnames = ["{}_marginal_vel_profile.png".format(x) for x in ('Ti', 'A', 'V', 'F')]
    factors = [1.0, 1.0, 1e-3, 1.0]

    if not use_plasma_F:
        xlabels.append('$\mathcal{F}$')
        ylabels.append('$P(\mathcal{F}) \Delta \mathcal{F}$')
        fnames.append("F_marginal_vel_profile.png")

    for idx, (xlabel, ylabel, factor) in enumerate(zip(xlabels, ylabels, factors)):
        print(np.mean(post[:,idx])*factor, np.std(post[:,idx])*factor)
        fig, ax = plt.subplots(figsize=(figsize, figsize/1.618))
        plotting.my_hist(ax, post[:, idx]*factor, bins='auto')
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        fig.savefig(path.join(plot_folder, fnames[idx]))
        # plt.close(fig)
        plt.show(block=True)
    

    for i, xlabel in enumerate(xlabels):
        for j in range(i+1, len(xlabels)):
            fig, ax = plt.subplots()
            cb = plotting.my_hist2d(ax, post[:, i]*factors[i], post[:, j]*factors[j])  
            cb.ax.tick_params(labelsize=fontsize)
            ax.set_xlabel(xlabels[j], fontsize=fontsize)
            ax.set_ylabel(xlabel, fontsize=fontsize)
            ax.tick_params(labelsize=fontsize)
            plt.show(block=False)
    plt.show()


def calculate_percentile_ranges(data):
    levels = [68, 95, 99]
    lower_levels = [50.0 - x/2.0 for x in levels]
    upper_levels = [50.0 + x/2.0 for x in levels]

    percentiles = {}
    for level, LL, UL in zip(levels, lower_levels, upper_levels):
        lower = np.percentile(data, LL, axis=0)
        upper = np.percentile(data, UL, axis=0)
        percentiles[level] = (lower, upper)

    return percentiles

def check_multi_solver(output_folder, locs, folders, Lpost, dpost, Fpost):

    figsize = (3.5, 3.5/1.618)
    fontsize = 9

    # locs = [5, 15, 25, 35]
    # folder = "/home/milhone/Research/python_FabryPerot/Data/PCX_Syn/"
    # folders = [path.join(folder, "{0:d}".format(x)) for x in locs]

    # output_folder = "/home/milhone/Research/python_FabryPerot/Data/PCX_Syn/multi"

    file_io.prep_folder(output_folder)
    plot_folder = path.join(output_folder, 'Plots')
    file_io.prep_folder(plot_folder)

    for folder in folders:
        print(folder)

    r_list = []
    s_list = []
    sd_list = []

    # read in data
    for folder in folders:
        data = file_io.h5_2_dict(path.join(folder, "argon_input.h5"))
        ix = data['fit_ix']['0']  # [0:-1:2]
        r_list.append(data['r'][ix])
        s_list.append(data['sig'][ix])
        sd_list.append(data['sig_sd'][ix])
    
    n_params = 7
    analyzer = pymultinest.Analyzer(n_params, outputfiles_basename=path.join(output_folder, "Ti_multi_"))
    post = analyzer.get_equal_weighted_posterior()
    nrows, ncols = post.shape

    # Let's sample from all of the posterior results
    n_samples = 200
    nL = len(Lpost)
    idx_L = np.random.choice(nL, size=n_samples)

    F_folder = "/home/milhone/Research/python_FabryPerot/Data/2018_10_28/"
    Fpost = np.loadtxt(path.join(F_folder, 'full_post_equal_weights.dat'), ndmin=2)[:,2]
    nF = len(Fpost)
    idx_F = np.random.choice(nF, size=n_samples)
    
    L_sampled = Lpost[idx_L]
    d_sampled = dpost[idx_L]
    #F_sampled = Fpost[idx_L]
    F_sampled = Fpost[idx_F]

    idx_cube = np.random.choice(nrows, size=n_samples)
    cubes = post[idx_cube, :]

    vals = []
    for r in r_list:
        vals.append(np.zeros((n_samples, len(r))))

    # Calculate the model values for all of the samples 
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # submit all of the calcultions to the executor and map to their index
        futures_map = dict()
        for j, (L, d, F, cube) in enumerate(zip(L_sampled, d_sampled, F_sampled, cubes)):
            Lnu = 100.0 * plasma.Lnu(cube[2], cube[0], mu=40, noise=False)
            fut = executor.submit(calculate_multi_models, r_list, locs, L, d, F, Lnu, cube)
            futures_map[fut] = j

        # grab the results from completed futures and store in vals
        for future in concurrent.futures.as_completed(futures_map):
            index = futures_map[future]
            output = future.result()
            for idx, out in enumerate(output):
                vals[idx][index, :] = out

    for idx, (r, sig, error) in enumerate(zip(r_list, s_list, sd_list)):
        fig, ax = plt.subplots(figsize=figsize)
        ax.fill_between(r, np.min(vals[idx], axis=0), np.max(vals[idx], axis=0), color='C0', alpha=0.5, label='Model')
        ax.errorbar(r, sig, yerr=error, errorevery=10, label='Data', color='C1')
        ax.set_title("Chord at r={0:3.1f} cm".format(locs[idx]))
        ax.legend(fontsize=fontsize)
        ax.set_xlabel("R (px)", fontsize=fontsize)
        ax.set_ylabel("Counts", fontsize=fontsize)
        fig.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.18)
        plt.show(block=False)
    plt.show(block=True)

    Lnu = np.zeros(nrows)
    rarr = np.linspace(0.0, 50.0, 1001.0)
    v_profile = np.zeros((nrows, len(rarr)))

    for idx, (nen0, Ti, Vouter) in enumerate(zip(post[:, 2], post[:, 0], post[:, 1])):
        Lnu[idx] = 100.0 * plasma.Lnu(nen0, Ti, mu=40, noise=False)
        v_profile[idx, :] = plasma.pcx_velocity_profile(rarr, Lnu[idx], 35.0, Vouter)

    fig, ax = plt.subplots(figsize=figsize)
    plotting.my_hist(ax, Lnu, bins='auto')
    print(np.mean(Lnu), np.std(Lnu))
    ax.axvline(24.816960868, color='k')
    ax.set_xlabel("$L_{\\nu}$ (cm)", fontsize=fontsize)
    ax.set_ylabel("$P(L_{\\nu}) \Delta L_{\\nu}$ (%)", fontsize=fontsize)
    fig.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.18)
    plt.show()

    Lnu_answer = 24.816960868
    r_outer = 35.0
    v_outer_answer = 1000.0
    vel_answer = plasma.pcx_velocity_profile(rarr, Lnu_answer, r_outer, v_outer_answer)

    v_perc = calculate_percentile_ranges(v_profile)
    levels = [68, 95, 99]
    alphas = [0.8, 0.5, 0.2]

    fig, ax = plt.subplots(figsize=figsize)
    for level, alpha in zip(levels, alphas):
        ax.fill_between(rarr, v_perc[level][0]/1000.0, v_perc[level][1]/1000.0, alpha=alpha, color='C0', label='{0:d}%'.format(level))
    mean_v = np.mean(v_profile, axis=0)
    #ax.fill_between(rarr, np.min(v_profile, axis=0), np.max(v_profile, axis=0), color='C0', alpha=0.5)
    #ax.plot(rarr, mean_v, 'C1', label='Mean')
    ax.plot(rarr, vel_answer/1000.0, color='C1', label='Answer', lw=1)
    ax.legend()
    ax.tick_params(labelsize=fontsize)
    ax.set_xlabel("R (cm)", fontsize=fontsize)
    ax.set_ylabel("V (km/s)", fontsize=fontsize)
    fig.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.18)
    fig.savefig(path.join(plot_folder, "vel_profile_prediction.pdf"), transparent=True)
    plt.show()

    # Make Marginal Posterior Plots
    xlabels = ['$T_i$ (eV)', '$V_{pk}$ (km/s)', '$n_e n_0$ ($10^{35}$ m${}^{-6}$)',
               '$A_0$ ($10^3$ Counts)', '$A_1$ ($10^3$ Counts)', '$A_2$ ($10^3$ Counts)',
               '$A_3$ ($10^3$ Counts)']

    ylabels = ['$P(T_i) \Delta T_i$ (%)', '$P(V_{pk}) \Delta V_{pk}$ (%)',
               '$P(n_e n_0) \Delta \left(n_e n_0\\right)$ (%)', '$P(A_0) \Delta A_0$ (%)',
               '$P(A_1) \Delta A_1$ (%)', '$P(A_2) \Delta A_2$ (%)', '$P(A_3) \Delta A_3$ (%)']

    fnames = ["{}_marginal_vel_profile.pdf".format(x) for x in ('Ti', 'V', 'nen0', 'A0', 'A1', 'A2', 'A3')]
    factors = [1.0, 1e-3, 1.e-35, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3]
    answers = [0.5, 1.0, 1.6, None, None, None, None]
    for idx, (xlabel, ylabel, factor) in enumerate(zip(xlabels, ylabels, factors)):
        fig, ax = plt.subplots(figsize=figsize)
        print(xlabel, np.mean(post[:, idx]*factor), np.std(post[:, idx]*factor))
        plotting.my_hist(ax, post[:, idx]*factor, bins='auto')
        ax.set_xlabel(xlabel, fontsize=fontsize)
        ax.set_ylabel(ylabel, fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        if answers[idx] is not None:
            ax.plot(answers[idx], 8.0, '*', ms=10, color='C3')
        fig.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.18)
        fig.savefig(path.join(plot_folder, fnames[idx]), transparent=True)
        # plt.close(fig)
        plt.show(block=True)

    for i, xlabel in enumerate(xlabels):
        for j in range(i+1, 3):
            fig, ax = plt.subplots()
            cb = plotting.my_hist2d(ax, post[:, i]*factors[i], post[:, j]*factors[j])  
            cb.ax.tick_params(labelsize=fontsize)
            ax.set_xlabel(xlabels[j], fontsize=fontsize)
            ax.set_ylabel(xlabel, fontsize=fontsize)
            ax.tick_params(labelsize=fontsize)
            plt.show(block=False)
    plt.show()



def calculate_multi_models(rlist, locs, L, d, F, Lnu, cube, nlambda=2000, nr=400):
    vals = []
    for idx, loc in enumerate(locs):
        w, spec = plasma.calculate_pcx_chord_emission(loc,
                         cube[0], w0, mu, Lnu, cube[1], nr=nr, nlambda=nlambda,
                         Lne=4.0, R_outer=35.0, rmax=42.0)

        vals.append(cube[idx+3] * models.general_model(rlist[idx], L, d, F, w, spec))

    return vals

