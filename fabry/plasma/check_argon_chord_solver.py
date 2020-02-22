from __future__ import print_function, division

import concurrent.futures
import os.path as path
import itertools

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pymultinest
from matplotlib import rcParams

import fabry.plasma.argon_chord_solver as acs
from fabry.plasma import plasma
from fabry.tools import plotting, file_io
import fabry.distributions.distribution as dist

rcParams['font.size'] = 16
rcParams['figure.figsize'] = (8, 6)

w0 = 487.98634
mu = 39.94

def prior(cube, ndim, nparams):
    """Transforms prior hypercube to physical units
    """
    # cube[0] = dist.UniformPrior(cube[0], 0.025, 3.0)
    # cube[1] = dist.UniformPrior(cube[1], -5000.0, 5000.0)
    cube[0] = dist.UniformPrior(cube[0], 0.025, 3.0)
    cube[1] = dist.UniformPrior(cube[1], 0.025, 3.0)
    cube[2] = dist.UniformPrior(cube[2], -2000.0, 2000.0)
    cube[3] = dist.UniformPrior(cube[3], -500, 500)

    # I want to use a less broad prior for the amplitudes, but not allow negative amplitude
    # for i in range(2, nparams):
    for i in range(3, nparams):
        #cube[i] = dist.TruncatedNormal(cube[i], 80.0, 30.0, 0.0, 1000.0) # 2020_02_05
        cube[i] = dist.TruncatedNormal(cube[i], 250.0, 150.0, 0.0, 1000.0) # 2020_01_15
        # cube[i] = dist.UniformPrior(cube[i], 1.0, 100.0)


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
    # vel_offsets = [0.0 for _ in range(n_images)]
    # vel_offsets = [322.0, 190.0, 162.0, 204.0]  # 2020_02_05
    vel_offsets = [322.0, 190.0, 162.0, 204.0, 306.0, 220.0]  # 2020_02_05 6 images
    # vel_offsets = [620.0, 470.0, 566.0, 509.0]  # 2020_01_15
    #vel_offsets = [620.0, 470.0, 566.0, 509.0, 503.0, 538.0]  # 2020_01_15  6 images

    # params are Ti, Vmax, +amplitudes for each image
    # n_params = 2 + n_images

    # params are Ti inner, Ti_outer, Vmax, +amplitudes for each image
    #n_params = 3 + n_images

    # params are Ti inner, Ti_outer, Vmax, Voffset, +amplitudes for each image
    n_params = 4 + n_images

    analyzer = pymultinest.Analyzer(n_params, outputfiles_basename=path.join(output_folder, "Ti_chord_"))
    trace = analyzer.get_data()[:,2::]
    prior(trace, None, n_params)

    post = analyzer.get_equal_weighted_posterior()
    post_means = np.mean(post, axis=0)
    post_std = np.std(post, axis=0)
    for i in range(n_params):
        print(f"{i}: {post_means[i]} +/- {post_std[i]}")
    nrows, ncols = post.shape

    alphas = [0.5, 0.3, 0.1]

    ###########################
    ## Velocity Profile Plot ##
    ###########################
    fig, ax = plt.subplots()
    linear_r = np.linspace(0.0, 40.0, 1000)
    # vprofile = plasma.pcx_velocity_profile(linear_r[None, :], 50.0, 32.0, post[:, 2][:, None]) 
    vprofile = plasma.pcx_velocity_profile(linear_r[None, :], 50.0, 32.0, post[:, 2][:, None], post[:, 3][:, None]) 
    levels = [68, 95, 99]
    vel_percentiles = calculate_percentile_ranges(vprofile, levels)
    vel_mean = np.mean(vprofile, axis=0)

    for alpha, level in zip(alphas, levels):
        ax.fill_between(linear_r, -vel_percentiles[level][0], 
                -vel_percentiles[level][1], color='k', alpha=alpha, 
                label=f"{level}%")
    ax.plot(linear_r, -vel_mean, 'C3', label='Mean')
    ax.legend(loc='upper left')
    ax.set_xlabel("R (cm)")
    ax.set_ylabel(r"$V_{\phi}$ (m/s)")
    fig.savefig(path.join(plot_folder, "V_profile.png"), dpi=400)
    plt.show(block=False)

    #####################
    ## Ti Profile Plot ##
    #####################
    Ti_profile = plasma.linear_Ti_profile(linear_r[None, :], post[:,0][:, None], post[:,1][:, None], 40.0) 
    Ti_percentiles = calculate_percentile_ranges(Ti_profile, levels)
    Ti_mean = np.mean(Ti_profile, axis=0)

    fig, ax = plt.subplots()
    for alpha, level in zip(alphas, levels):
        ax.fill_between(linear_r, Ti_percentiles[level][0], 
                Ti_percentiles[level][1], color='k', alpha=alpha, 
                label=f"{level}%")
    ax.plot(linear_r, Ti_mean, 'C3', label='Mean')
    ax.legend(loc='lower right')
    ax.set_xlabel("R (cm)")
    ax.set_ylabel(r"$T_{i}$ (eV)")
    _, yhigh = ax.get_ylim()
    ax.set_ylim(-.02, 1.15*yhigh)
    fig.savefig(path.join(plot_folder, "Ti_profile.png"), dpi=400)
    plt.show()


    idx_cube = np.random.choice(nrows, size=n_samples)
    cubes = post[idx_cube, :]

    # fig, ax = plt.subplots()
    # sns.kdeplot(post[:, 0], post[:, 1], shade=True, shade_lowest=False, cmap='Reds')
    # plt.show()

    # allocate memory for posterior values
    values = []
    for rr in r:
        values.append(np.zeros((n_samples, len(rr))))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures_map = dict()
        # for j, (L, d, F, cube) in enumerate(zip(L_sampled, d_sampled, F_sampled, cubes)):
        for j, cube in enumerate(cubes):
            fut = executor.submit(calculate_multi_models, r, impact_factors, vel_offsets, L, d, F, Lnu, cube)
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
    # for i, (rr, ss, std, val_percentile) in enumerate(zip(r, sig, sd, val_percentiles)):
    #     fig, ax = plt.subplots()
    #     for alpha, level in zip(alphas, levels):
    #         ax.fill_between(rr, val_percentile[level][0], val_percentile[level][1], color='k', alpha=alpha)
    #     #ax.fill_between(rr, values[i].min(axis=0), values[i].max(axis=0), color='k', alpha=0.3)
    #     ax.plot(rr, val_mean[i], color='k', label='Fit')
    #     # ax.errorbar(rr, ss, yerr=std, color='C1', errorevery=1, label='Data')
    #     ax.fill_between(rr, ss-std, ss+std, color='C1', alpha=0.5)
    #     ax.plot(rr, ss, color='C1', label='Data')
    #     impact_factor = impact_factors[i]

    #     title_string = f"R = {impact_factor:2.1f} cm"
    #     ax.set_title(title_string)
    #     ax.set_xlabel("R (px)")
    #     ax.set_ylabel("Counts")
    #     ax.legend()
    #     plt.show(block=False)
    # plt.show(block=True)

    ############################
    ## Fit Plot with Subplots ##
    ############################

    # Figure out grid shape
    n_images = len(r)
    nrows = 2
    ncols, remainder = divmod(n_images, nrows)
    if remainder != 0:
        nrows += 1
    alphas = [0.5, 0.3, 0.1]

    # Let's reorder the plots by impact factor
    print(impact_factors)
    isort = np.argsort(impact_factors)
    print(isort)
    impact_factors = [impact_factors[x] for x in isort]
    print(impact_factors)
    r = [r[x] for x in isort]
    sig = [sig[x] for x in isort]
    sd = [sd[x] for x in isort]
    val_percentiles = [val_percentiles[x] for x in isort]
    val_mean = [val_mean[x] for x in isort]

    fig, ax = plt.subplots(nrows, ncols, sharex=True)
    for i, (rr, ss, std, val_percentile) in enumerate(zip(r, sig, sd, val_percentiles)):
        r, c = divmod(i, ncols)

        # plot fit with percentiles
        for alpha, level in zip(alphas, levels):
            ax[r,c].fill_between(rr, val_percentile[level][0], val_percentile[level][1], color='k', alpha=alpha)
        ax[r,c].plot(rr, val_mean[i], color='k', label='Fit')

        # plot data with error (via fill between)
        print(rr.shape, ss.shape, std.shape)
        ax[r,c].fill_between(rr, ss-std, ss+std, color='C1', alpha=0.5)
        ax[r,c].plot(rr, ss, color='C1', label='Data')

        impact_factor = impact_factors[i]
        ylow, yhigh = ax[r,c].get_ylim()
        ax[r,c].set_ylim([ylow, 1.25*yhigh])
        xleft, _ = ax[r,c].get_xlim()
        text_string = f"{impact_factor:2.1f} cm"
        xy_loc = (1.01*xleft, 1.1*yhigh)
        ax[r,c].annotate(text_string, xy_loc, ha='left', va='top')

    # legend placement logic
    if ncols %2 == 0:
        # even
        i = ncols // 2 - 1
        ax[0,i].legend(ncol=2, loc='upper center', bbox_to_anchor=(1.0+0.1, 1.25))
    else:
        i = ncols // 2 
        ax[0,i].legend(ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.25))

    # axis labels
    for i in range(ncols):
        ax[nrows-1, i].set_xlabel("R (px)")
    for i in range(nrows):
        ax[i, 0].set_ylabel("Counts")

    # plot tidyness
    fig.tight_layout()
    fig.subplots_adjust(hspace=0.0, wspace=0.2)
    for aa in ax:
        for a in aa:
            low, high = a.get_ylim()
            yticks = a.get_yticks()
            a.set_yticks(yticks[0:-1])
            a.set_ylim(low, high)
    fig.savefig(path.join(plot_folder, "ringsum_fits_all_chords.png"), dpi=400)
    plt.show()

    ######################## 
    ## Histogram Plotting ##
    ######################## 
    #xlabels = [r'$T_{i,\rm{center}}$ (eV)', r'$T_{i,\rm{edge}}$ (eV)', '$V_{pk}$ (km/s)']
    xlabels = [r'$T_{i,\rm{center}}$ (eV)', r'$T_{i,\rm{edge}}$ (eV)', '$V_{pk}$ (km/s)', r"$V_{\rm{offset}}$ (km/s)"]
    #ylabels = [r'PDF$(T_{i,\rm{center}})$', r'PDF$(T_{i,\rm{edge}})$', 'PDF$(V_{pk})$']
    ylabels = [r'PDF$(T_{i,\rm{center}})$', r'PDF$(T_{i,\rm{edge}})$', 'PDF$(V_{pk})$', r"PDF$(V_{\rm{offset}})$"]
    factors = [1.0, 1.0, 1e-3, 1e-3]
    #post_variables = ["Ti_center", "Ti_edge", "Vpeak",]
    post_variables = ["Ti_center", "Ti_edge", "Vpeak", "Voffset"]
    post_variables += [f"A{i}" for i in range(n_images)]
    marginal_fnames = [path.join(plot_folder, f"{variable}_marginal.png") for variable in post_variables]
    # add in the amplitudes
    xlabels += ["$A_{0}$".format(i) for i in range(n_images)]
    ylabels += ["PDF$(A_{0})$".format(i) for i in range(n_images)]
    factors += [1.0 for _ in range(n_images)]

    for i in range(n_params):
        fig, ax = plt.subplots()
        sns.distplot(post[:, i]*factors[i], kde=True, color='C3')
        ax.set_xlabel(xlabels[i])
        ax.set_ylabel(ylabels[i])
        fig.savefig(marginal_fnames[i], dpi=400)
        plt.show(block=False)

    plt.show(block=True)

    #make_histograms(post, xlabels, ylabels, factors, None, save=False, block=False)
    #plt.show()
    ##################################################### 
    ## 2D Marginal Posterior Distribution for Ti and V ##
    ##################################################### 
    #fig, ax = plt.subplots()
    #sns.kdeplot(post[:, 1]*factors[1], post[:, 0]*factors[0], cmap='Reds', shade=True, shade_lowest=False)
    ##cb = plotting.my_hist2d(ax, post[:, 0]*factors[0], post[:, 1]*factors[1], plot_contours=True)
    ##cb.ax.tick_params(labelsize=16)
    #ax.set_xlabel(xlabels[1])
    #ax.set_ylabel(xlabels[0])
    #plt.show(block=False)

    #fig, ax = plt.subplots()
    #sns.kdeplot(post[:, 2]*factors[2], post[:, 0]*factors[0], cmap='Reds', shade=True, shade_lowest=False)
    ##cb = plotting.my_hist2d(ax, post[:, 0]*factors[0], post[:, 2]*factors[2], plot_contours=True)
    ##cb.ax.tick_params(labelsize=16)
    #ax.set_xlabel(xlabels[2])
    #ax.set_ylabel(xlabels[0])
    #plt.show(block=False)

    #fig, ax = plt.subplots()
    #sns.kdeplot(post[:, 2]*factors[2], post[:, 1]*factors[1], cmap='Reds', shade=True, shade_lowest=False)
    ##cb = plotting.my_hist2d(ax, post[:, 1]*factors[1], post[:, 2]*factors[2], plot_contours=True)
    ##cb.ax.tick_params(labelsize=16)
    #ax.set_xlabel(xlabels[2])
    #ax.set_ylabel(xlabels[1])
    #plt.show()


    # fig, ax = plt.subplots(2, 2)
    # ax[0, 1].axis('off')
    # sns.distplot(post[:, 1], kde=True, color='C3', ax=ax[0, 0])
    # sns.distplot(post[:, 0], kde=True, color='C3', ax=ax[1, 1], vertical=True)
    # sns.kdeplot(post[:, 1]*factors[1], post[:, 0]*factors[0], cmap='Reds', shade=True, shade_lowest=False, ax=ax[1, 0])

    # ax[0, 0].set_xticks([])
    # ax[1, 1].set_yticks([])
    # ax[1, 0].set_xlabel(xlabels[1])
    # ax[1, 0].set_xlabel(xlabels[0])

    # fig.subplots_adjust(hspace=0.0, wspace=0.0)
    # plt.show(block=False)

    #g = sns.jointplot(post[:, 0], post[:, 1], kind='kde', color='C3', join_kws={'shade_lowest': False})

    n = 4
    pairs = itertools.product(range(n), range(n))
    pairs = filter(lambda tup: tup[1] > tup[0], pairs)
    for (i, j) in pairs:
        g = sns.JointGrid(post[:, i], post[:, j], space=0.0, height=8)
        g = g.plot_joint(sns.kdeplot, cmap='Reds', shade_lowest=False, shade=True)
        g = g.plot_marginals(sns.distplot, kde=True, color='C3')
        g.set_axis_labels(xlabels[i], xlabels[j])
        fig = g.fig
        fig.tight_layout()
        fig.subplots_adjust(hspace=0.0, wspace=0.0)
        fname = f"{post_variables[i]}_{post_variables[j]}_joint_marginal.png"
        fname = path.join(plot_folder, fname)
        fig.savefig(fname, dpi=400)
        plt.show(block=False)

    # g = sns.JointGrid(post[:, 0], post[:, 1], space=0.0, height=8)
    # g = g.plot_joint(sns.kdeplot, cmap='Reds', shade_lowest=False, shade=True)
    # g = g.plot_marginals(sns.distplot, kde=True, color='C3') 
    # g.set_axis_labels(xlabels[0], xlabels[1])
    # fig = g.fig
    # fig.tight_layout()
    # fig.subplots_adjust(hspace=0.0, wspace=0.0)
    # plt.show(block=False)

    # g = sns.JointGrid(post[:, 0], post[:, 2]*factors[2], space=0.0, height=8)
    # g = g.plot_joint(sns.kdeplot, cmap='Reds', shade_lowest=False, shade=True)
    # g = g.plot_marginals(sns.distplot, kde=True, color='C3') 
    # g.set_axis_labels(xlabels[0], xlabels[2])
    # fig = g.fig
    # fig.tight_layout()
    # fig.subplots_adjust(hspace=0.0, wspace=0.0)
    # plt.show(block=False)

    # g = sns.JointGrid(post[:, 1], post[:, 2]*factors[2], space=0.0, height=8)
    # g = g.plot_joint(sns.kdeplot, cmap='Reds', shade_lowest=False, shade=True)
    # g = g.plot_marginals(sns.distplot, kde=True, color='C3') 
    # g.set_axis_labels(xlabels[1], xlabels[2])
    # fig = g.fig
    # fig.tight_layout()
    # fig.subplots_adjust(hspace=0.0, wspace=0.0)
    # plt.show(block=False)

    plt.show()


def make_histograms(posterior, xlabels, ylabels, factors, fnames, save=False, plot_folder=None, block=True):

    for idx, (xlabel, ylabel, factor) in enumerate(zip(xlabels, ylabels, factors)):
        hist_data = posterior[:, idx]*factor
        fig, ax = plt.subplots()
        bins = plotting.my_hist(ax, hist_data, bins='auto')

        mu = np.mean(hist_data)
        sigma = np.std(hist_data)
        x = np.linspace(hist_data.min(), hist_data.max(), 500) 
        kde = np.exp(-(x-mu)**2/2/sigma**2) / np.sqrt(2.0 *np.pi) / sigma
        ax.plot(x, kde*(bins[1]-bins[0])*100, zorder=100, color='k')
        # print("")
        # print(xlabel)
        # print('bin width', bins[1]-bins[0])
        # print('kde max', kde.max())
        # print("")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        if save and plot_folder is not None:
            filepath = path.join(plot_folder, fnames[idx])
            fig.savefig(filepath, transparent=True)

        plt.show(block=block)

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
        V_params = [cube[2], Lnu, r_anode, cube[3]]
        ne_params = [r_anode, ]
        delta_d = acs.calculate_delta_d(v_offset)
        # amp = cube[idx + 3]
        amp = cube[idx + 4]

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








