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
                        resume=resume, verbose=True, sampling_efficiency='model', n_live_points=200,
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

    def log_likelihood(cube, ndim, nparams):
        iL = np.random.choice(nL)
        LL = Lpost[iL]
        dd = dpost[iL]
        jF = np.random.choice(nF)
        FF = Fpost[jF]

        #print(LL, dd, FF)
        vals = models.forward_model(r, LL, dd, FF, w0, mu, cube[1], cube[0], cube[2], nlambda=2000)

        chisq = np.sum((vals - sig)**2 / error**2)
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

    v_lim = [-10000, 10000]

    nL = len(Lpost)
    nF = len(Fpost)
    n_params = 3

    if test_plot:
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
                        resume=resume, verbose=True, sampling_efficiency='model', n_live_points=200,
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
        cube[3] = 10 ** (cube[3] * (Lnu_lim[1] - Lnu_lim[0]) + Lnu_lim[0])

    def log_likelihood(cube, ndim, nparams):
        iL = np.random.choice(nL)
        jF = np.random.choice(nF)
        LL = Lpost[iL]
        dd = dpost[iL]
        FF = Fpost[jF]

        w, spec = plasma.calculate_pcx_chord_emission(impact_factor,
                                                      cube[0], w0, mu, cube[3], cube[2], nr=nr, nlambda=nlambda, Lne=4.0, R_outer=35.0, rmax=42.0)

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

    Lnu_lim = [0.1, 100]
    Lnu_lim = [np.log10(lim) for lim in Lnu_lim]

    nL = len(Lpost)
    nF = len(Fpost)

    n_params = 4

    impact_factor = 33.0
    nr = 400
    nlambda = 2000

    if test_plot:
        # do a test plot
        npts = 100
        test_sig = np.zeros((npts, len(r)))
        output = np.loadtxt(path.join(output_folder, 'Ti_profileV_post_equal_weights.dat'), ndmin=2)[:, 0:-1]
        nout = output.shape[0]
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

            w, spec = plasma.calculate_pcx_chord_emission(impact_factor,
                                                          cub[0], w0, mu, cub[3], cub[2], nr=nr, nlambda=nlambda, R_outer=35.0, Lne=4.0, rmax=42.0)

            test_sig[idx, :] = cub[1] * models.general_model(r, L, d, F, w, spec)

        fig, ax = plt.subplots()
        for i in range(npts):
            ax.plot(r, test_sig[i, :], 'C0')
        ax.errorbar(r, sig, yerr=error, color='C1')
        w, pcx_emis = plasma.calculate_pcx_chord_emission(impact_factor, 0.5, w0, mu, 20.0, 4000.0, rmax=42.0, nr=500, nlambda=2000, Lne=4.0, R_outer=35.0)
        vals = 1.01*np.mean(output[:, 1]) * models.general_model(r, 150.0/0.004, 0.88, 20.7, w, pcx_emis)
        #ax.plot(r, vals, 'C2')
        plt.show()

    else:
        # run multinest
        pymultinest.run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
                        resume=resume, verbose=True, sampling_efficiency='model', n_live_points=500,
                        outputfiles_basename=path.join(output_folder, 'Ti_profileV_'))
