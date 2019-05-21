from __future__ import division, print_function
from ..core import models
from ..tools import file_io
import plasma
import os.path as path
import numpy as np
import pymultinest
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass


w0 = 487.98634
mu = 39.948


def no_vel_solver(output_folder, Fpost, Lpost, dpost, resume=True, test_plot=False):
    """PyMultinest Solver for Ar II with no velocity

    Args:
        output_folder (str): path to folder for input and output folders
        Fpost (np.ndarray): posterior results for the finesse
        Lpost (np.ndarray): posterior results for the camera focal length
        dpost (np.ndarray): posterior results for the etalon spacing
        resume (bool): resume calculation if True, default=True
        test_plot (bool): make a test plot of prior intstead of solving, default=False
    """
    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0] * (Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]
        cube[1] = cube[1] * (A_lim[1] - A_lim[0]) + A_lim[0]

    def log_likelihood(cube, ndim, nparams):
        iL = np.random.choice(nL)
        # jF = np.random.choice(nF)
        L = Lpost[iL]
        d = dpost[iL]
        # F = Fpost[jF]
        F = Fpost[iL]

        vals = models.forward_model(r, L, d, F, w0, mu, cube[1], cube[0], 0.0, nlambda=2000)

        chisq = np.sum((vals - sig) ** 2 / error ** 2)
        return -chisq / 2

    data_filename = path.join(output_folder, "argon_input.h5")
    data = file_io.h5_2_dict(data_filename)

    ix = data['fit_ix']['0'][0:-1:2]
    r = data['r'][ix]
    sig = data['sig'][ix]
    error = data['sig_sd'][ix]

    A_max = np.max(sig)
    A_lim = [0.5 * A_max, 2.0 * A_max]

    Ti_lim = [0.025, 3.0]

    nL = len(Lpost)
    nF = len(Fpost)
    n_params = 2

    if test_plot:
        # do a test plot
        npts = 100
        test_sig = np.zeros((npts, len(r)))
        for idx in range(npts):
            cub = np.random.random(size=n_params)
            log_prior(cub, None, None)
            i = np.random.choice(nL)
            j = np.random.choice(nF)
            LL = Lpost[i]
            dd = dpost[i]
            FF = Fpost[j]
            test_sig[idx, :] = models.forward_model(r, LL, dd, FF, w0, mu, cub[1], cub[0],
                                                  0.0, nlambda=2000)

        fig, ax = plt.subplots()
        for i in range(npts):
            ax.plot(r, test_sig[i, :], 'C0')
        ax.errorbar(r, sig, yerr=error, color='C1')
        plt.show()

    else:
        # run multinest
        pymultinest.run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
                        resume=resume, verbose=True, sampling_efficiency='model', n_live_points=400,
                        outputfiles_basename=path.join(output_folder, 'Ti_noV_'))


def const_vel_solver(output_folder, Fpost, Lpost, dpost, resume=True, test_plot=False):
    """PyMultinest Solver for Ar II with constant velocity

    Args:
        output_folder (str): path to folder for input and output folders
        Fpost (np.ndarray): posterior results for the finesse
        Lpost (np.ndarray): posterior results for the camera focal length
        dpost (np.ndarray): posterior results for the etalon spacing
        resume (bool): resume calculation if True, default=True
        test_plot (bool): make a test plot of prior intstead of solving, default=False
    """
    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0] * (Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]
        cube[1] = cube[1] * (A_lim[1] - A_lim[0]) + A_lim[0]
        cube[2] = cube[2] * (v_lim[1] - v_lim[0]) + v_lim[0]
        #cube[3] = cube[3] * (F_lim[1] - F_lim[0]) + F_lim[0]
        #cube[3] = cube[3] * (off_lim[1] - off_lim[0]) + off_lim[0]
        #cube[4] = cube[4] * (F_lim[1] - F_lim[0]) + F_lim[0]

    def log_likelihood(cube, ndim, nparams):
        iL = np.random.choice(nL)
        LL = Lpost[iL]
        dd = dpost[iL]
        FF = Fpost[iL]
        #jF = np.random.choice(nF)
        #FF = Fpost[jF]
        #LL = 0.124073186132456485E+05
        #dd = 0.876551993778090566E+00
        #FF = 0.205637295793309107E+02
        #print(LL, dd, FF)
        #vals = models.forward_model(r, LL, dd, cube[3], w0, mu, cube[1], cube[0], cube[2], nlambda=2000)
        vals = models.forward_model(r, LL, dd, FF, w0, mu, cube[1], cube[0], cube[2], nlambda=2000)
        #vals = models.forward_model(r, LL, dd, FF, w0, mu, cube[1], cube[0], cube[2], nlambda=2000)
        #vals += cube[3]
        #vals += 2.7
        chisq = np.sum((vals - sig)**2 / error**2)
        return -chisq / 2

    #print('scaling L from npix=3 to npix=1')
    #Lpost = 3.0 * Lpost

    data_filename = path.join(output_folder, "argon_input.h5")
    data = file_io.h5_2_dict(data_filename)

    ix = data['fit_ix']['0']#[0:-1:2]
    r = data['r'][ix]
    sig = data['sig'][ix]
    error = data['sig_sd'][ix]

    A_max = np.max(sig)
    A_lim = [0.5 * A_max, 2.0 * A_max]

    # Ti_lim = [0.025, 3.0]
    Ti_lim = [0.0, 3.0]

    v_lim = [-10000, 10000]
    F_lim = [18, 35]
    #F_folder = "/home/milhone/Research/python_FabryPerot/Data/PCX/2158/"
    #Fpost = np.loadtxt(path.join(F_folder, 'Ti_constV_post_equal_weights.dat'), ndmin=2)[:,3]
    #F_folder = "/home/milhone/Research/python_FabryPerot/Data/2018_10_28/"
    #Fpost = np.loadtxt(path.join(F_folder, 'full_post_equal_weights.dat'), ndmin=2)[:,2]
    off_lim = [-20.0, 20.0]
    nL = len(Lpost)
    nF = len(Fpost)
    n_params = 3

    if False:#test_plot:
        # do a test plot
        npts = 100
        test_sig = np.zeros((npts, len(r)))
        print(test_sig.shape, r.shape)
        output = np.loadtxt(path.join(output_folder, 'Ti_constV_post_equal_weights.dat'), ndmin=2)[:, 0:-1]
        nout = output.shape[0]
        for idx in range(npts):
            k = np.random.choice(nout)
            #cub = np.random.random(size=n_params)
            cub = output[k, :]
            #log_prior(cub, None, None)
            i = np.random.choice(nL)
            j = np.random.choice(nF)
            L = Lpost[i]
            d = dpost[i]
            F = Fpost[j]
            test_sig[idx, :] = models.forward_model(r, L, d, F, w0, mu, cub[1], cub[0],
                                                  cub[2], nlambda=2000)
            #test = log_likelihood(cub, None, None)

        fig, ax = plt.subplots()
        for i in range(npts):
            ax.plot(r, test_sig[i, :], 'C0')
        ax.errorbar(r, sig, yerr=error, color='C1')
        plt.show()

    else:
        # run multinest
        pymultinest.run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
                        resume=resume, verbose=True, sampling_efficiency='model', n_live_points=400,
                        outputfiles_basename=path.join(output_folder, 'Ti_constV_'))

def profile_vel_solver(output_folder, Fpost, Lpost, dpost, resume=True, test_plot=False):
    """PyMultinest Solver for Ar II with velocity profile for PCX-U with outer boundary
    spinning

    Args:
        output_folder (str): path to folder for input and output folders
        Fpost (np.ndarray): posterior results for the finesse
        Lpost (np.ndarray): posterior results for the camera focal length
        dpost (np.ndarray): posterior results for the etalon spacing
        resume (bool): resume calculation if True, default=True
        test_plot (bool): make a test plot of prior intstead of solving, default=False
    """
    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0] * (Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]
        cube[1] = 10 ** (cube[1] * (A_lim[1] - A_lim[0]) + A_lim[0])
        cube[2] = cube[2] * (v_lim[1] - v_lim[0]) + v_lim[0]
        #cube[3] = 10 ** (cube[3] * (Lnu_lim[1] - Lnu_lim[0]) + Lnu_lim[0])
        cube[3] = 10 ** (cube[3] * (nen0_loglim[1] - nen0_loglim[0]) + nen0_loglim[0])

    def log_likelihood(cube, ndim, nparams):
        iL = np.random.choice(nL)
        jF = np.random.choice(nF)
        LL = Lpost[iL]
        dd = dpost[iL]
        FF = Fpost[jF]
        Lnu = 100.0 * plasma.Lnu(cube[3], cube[0], mu=40, noise=False)
        # print(Lnu)
        w, spec = plasma.calculate_pcx_chord_emission(impact_factor,
                                                      cube[0], w0, mu, Lnu, cube[2], nr=nr, nlambda=nlambda, Lne=4.0, R_outer=35.0, rmax=42.0)

        vals = cube[1] * models.general_model(r, LL, dd, FF, w, spec)

        chisq = np.sum((vals - sig) ** 2 / error ** 2)
        return -chisq / 2

    data_filename = path.join(output_folder, "argon_input.h5")
    data = file_io.h5_2_dict(data_filename)

    ix = data['fit_ix']['0'][0:-1:2]
    r = data['r'][ix]
    sig = data['sig'][ix]
    error = data['sig_sd'][ix]

    # I don't understand my integration for the emission very well
    # so I'm going to a huge prior in log space
    A_max = np.max(sig)
    A_lim = [0.1 * A_max, 10.0 * A_max]
    A_lim = [np.log10(AA) for AA in A_lim]

    Ti_lim = [0.025, 3.0]

    v_lim = [-10000, 10000]

    # Lnu_lim = [0.1, 100]
    # Lnu_lim = [np.log10(lim) for lim in Lnu_lim]

    # ne = [1e16, 1e19], n0 = [1e16, 1e19]
    # nen0_loglim = [16+16, 19+19]
    nen0_loglim = [np.log10(1e17*1e17), np.log10(2e18*2e18)]

    nL = len(Lpost)
    nF = len(Fpost)

    n_params = 4

    # impact_factor = 33.0
    impact_factor = 25.0
    nr = 400
    nlambda = 2000

    if test_plot:
        # do a test plot
        npts = 5000
        test_sig = np.zeros((npts, len(r)))
        # output = np.loadtxt(path.join(output_folder, 'Ti_profileV_post_equal_weights.dat'), ndmin=2)[:, 0:-1]
        Lnu_vals = np.zeros(npts)
        # nout = output.shape[0]
        for idx in range(npts):
            cub = np.random.random(size=n_params)
            log_prior(cub, None, None)
            #k = np.random.choice(nout)
            #cub = output[k, :]
            i = np.random.choice(nL)
            j = np.random.choice(nF)
            L = Lpost[i]
            d = dpost[i]
            F = Fpost[j]

            Lnu = 100.0 * plasma.Lnu(cub[3], cub[0], mu=40, noise=True)
            Lnu_vals[idx] = Lnu
            #w, spec = plasma.calculate_pcx_chord_emission(impact_factor,
            #                                              cub[0], w0, mu, Lnu, cub[2], nr=nr, nlambda=nlambda, R_outer=35.0, Lne=4.0, rmax=42.0)

            #test_sig[idx, :] = cub[1] * models.general_model(r, L, d, F, w, spec)

        #fig, ax = plt.subplots()
        #for i in range(npts):
        #    ax.plot(r, test_sig[i, :], 'C0')
        #ax.errorbar(r, sig, yerr=error, color='C1')
        ##w, pcx_emis = plasma.calculate_pcx_chord_emission(impact_factor, 0.5, w0, mu, 20.0, 4000.0, rmax=42.0, nr=500, nlambda=2000, Lne=4.0, R_outer=35.0)
        ## vals = 1.01*np.mean(output[:, 1]) * models.general_model(r, 150.0/0.004, 0.88, 20.7, w, pcx_emis)
        ##ax.plot(r, vals, 'C2')
        #plt.show()
        fig, ax = plt.subplots()
        ax.hist(Lnu_vals)
        plt.show()

    else:
        # run multinest
        pymultinest.run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
                        resume=resume, verbose=True, sampling_efficiency='model', n_live_points=500,
                        outputfiles_basename=path.join(output_folder, 'Ti_profileV_'))


def multi_image_solver(output_folder, locs, folders, Lpost, dpost, Fpost, test_plot=False, resume=True):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0] * (Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]
        cube[1] = cube[1] * (v_lim[1] - v_lim[0]) + v_lim[0]
        #cube[2] = 10 ** (cube[2] * (nen0_loglim[1] - nen0_loglim[0]) + nen0_loglim[0])
        cube[2] = cube[2] * (Lnu_lim[1] - Lnu_lim[0]) + Lnu_lim[0]
        cube[3] = 10 ** (cube[3] * (A_lim[0][1] - A_lim[0][0]) + A_lim[0][0])
        cube[4] = 10 ** (cube[4] * (A_lim[1][1] - A_lim[1][0]) + A_lim[1][0])
        cube[5] = 10 ** (cube[5] * (A_lim[2][1] - A_lim[2][0]) + A_lim[2][0])
        cube[6] = 10 ** (cube[6] * (A_lim[3][1] - A_lim[3][0]) + A_lim[3][0])

    def log_likelihood(cube, ndim, nparams):
        # pick L, d and F
        iL = np.random.choice(nL)
        jF = np.random.choice(nF)
        LL = Lpost[iL]
        dd = dpost[iL]
        FF = Fpost[jF]

        # calculate Lnu from Ti and nen0
        #Lnu = 100.0 * plasma.Lnu(cube[2], cube[0], mu=40, noise=False)

        # print(Lnu)

        # loop over the different chord locations and calculate chi squared
        chisq = 0.0
        for idx, (r, sig, error) in enumerate(zip(r_list, s_list, sd_list)):
            w, spec = plasma.calculate_pcx_chord_emission(locs[idx],
                               #cube[0], w0, mu, Lnu, cube[1], nr=nr, nlambda=nlambda, Lne=4.0, R_outer=35.0, rmax=42.0)
                               cube[0], w0, mu, cube[2], cube[1], nr=nr, nlambda=nlambda, Lne=4.0, R_outer=35.0, rmax=42.0)

            vals = cube[idx+3] * models.general_model(r, LL, dd, FF, w, spec)

            chisq += np.sum((vals - sig) ** 2 / error ** 2)

        return -chisq / 2.0


    # locs = [5, 15, 25, 35]
    # folder = "/home/milhone/Research/python_FabryPerot/Data/PCX_Syn/"
    # folders = [path.join(folder, "{0:d}".format(x)) for x in locs]

    # output_folder = "/home/milhone/Research/python_FabryPerot/Data/PCX_Syn/multi"
    file_io.prep_folder(output_folder)

    for folder in folders:
        print(folder)

    r_list = []
    s_list = []
    sd_list = []

    # read in data
    for folder in folders:
        data = file_io.h5_2_dict(path.join(folder, "argon_input.h5"))
        ix = data['fit_ix']['0'][0:-1:2]
        r_list.append(data['r'][ix])
        s_list.append(data['sig'][ix]-0.8)
        sd_list.append(data['sig_sd'][ix])

    # make amplitude prior limits
    A_lim = []
    for sig in s_list:
        A_max = np.max(sig)
        temp = [0.1 * A_max, 10.0 * A_max]
        A_lim.append([np.log10(AA) for AA in temp])

    Ti_lim = [0.025, 3.0]

    v_lim = [-10000, 10000]

    Lnu_lim = [1., 50.]
    # Lnu_lim = [np.log10(lim) for lim in Lnu_lim]

    # ne = [1e16, 1e19], n0 = [1e16, 1e19]
    # nen0_loglim = [16+16, 19+19]
    nen0_loglim = [np.log10(1e17*1e17), np.log10(2e18*2e18)]

    F_folder = "/home/milhone/Research/python_FabryPerot/Data/2018_10_28/"
    Fpost = np.loadtxt(path.join(F_folder, 'full_post_equal_weights.dat'), ndmin=2)[:,2]
    Lpost *= 3.0
    nL = len(Lpost)
    nF = len(Fpost)

    n_params = 7

    nr = 400
    nlambda = 2000

    if test_plot:
        pass
    else:
        # run multinest
        pymultinest.run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
                        resume=resume, verbose=True, sampling_efficiency='model', n_live_points=600,
                        outputfiles_basename=path.join(output_folder, 'Ti_multi_'))


