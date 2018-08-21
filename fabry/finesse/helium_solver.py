from __future__ import division, print_function
import os.path as path
from ..tools import file_io as io
import matplotlib.pyplot as plt
import numpy as np
from ..core import models
from pymultinest import run

w0 = 468.619458
mu = 232.03806


def full_solver(output_folder, data_filename, resume=True, test_plot=False):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(L_lim[1] - L_lim[0]) + L_lim[0]
        cube[1] = cube[1]*(d_lim[1] - d_lim[0]) + d_lim[0]
        cube[2] = cube[2]*(F_lim[1] - F_lim[0]) + F_lim[0]
        cube[3] = cube[3]*(A_lim[1] - A_lim[0]) + A_lim[0]
        cube[4] = cube[4]*(B_lim[1] - B_lim[0]) + B_lim[0]
        cube[5] = cube[5]*(Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]

    def log_likelihood(cube, ndim, nparams):
        vals0, vals1 = forward_model(cube)
        chisq = np.sum((vals0 - sig0)**2 / sig0_sd**2)
        chisq += np.sum((vals1 - sig1)**2 / sig1_sd**2)

        return -chisq / 2.0

    def forward_model(cube):
        vals0 = models.offset_forward_model(r0, cube[0], cube[1], cube[2], w0,
                                            mu, cube[3], cube[5], 0.0, coeff=0.05)
        vals1 = models.offset_forward_model(r1, cube[0], cube[1], cube[2], w0,
                                            mu, cube[4], cube[5], 0.0, coeff=0.05)
        return vals0, vals1

    data = io.h5_2_dict(data_filename)
    ix0 = data['fit_ix']['0']
    ix1 = data['fit_ix']['1']

    r0 = data['r'][ix0]
    sig0 = data['sig'][ix0]
    sig0_sd = data['sig_sd'][ix0]

    r1 = data['r'][ix1]
    sig1 = data['sig'][ix1]
    sig1_sd = data['sig_sd'][ix1]

    L_lim = [135.0, 150.0]
    L_lim = [x / 0.004 for x in L_lim]

    d_lim = [0.7, 0.9]

    F_lim = [17.0, 21.0]

    A_max = np.max(sig0)
    A_lim = [0.75*A_max, 2.0*A_max]

    B_max = np.max(sig1)
    B_lim = [0.75*B_max, 2.0*B_max]

    Ti_lim = [0.025, 0.3]

    n_params = 6
    folder = path.abspath(output_folder)

    if test_plot:
        print('*****************')
        print('*   Test Plot   *')
        print('*****************')
        npts = 100
        test0 = np.zeros((npts, len(r0)))
        test1 = np.zeros((npts, len(r1)))
        for i in range(npts):
            cube = [np.random.random() for _ in range(n_params)]
            log_prior(cube, None, None)
            test0[i, :], test1[i, :] = forward_model(cube)

        fig, ax = plt.subplots()
        ax.errorbar(r0, sig0, yerr=sig0_sd, color='C1')
        ax.errorbar(r1, sig1, yerr=sig1_sd, color='C1')
        for i in range(npts):
            ax.plot(r0, test0[i, :], 'C0')
            ax.plot(r1, test1[i, :], 'C0')
        plt.show()

    else:
        run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
            resume=resume, verbose=True, sampling_efficiency='q', n_live_points=75,
            outputfiles_basename=path.join(folder, 'full_'))


def solver(output_folder, prior_filename, data_filename, Lpost, dpost, resume=True, test_plot=True):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(F_lim[1] - F_lim[0]) + F_lim[0]
        cube[1] = cube[1]*(A_lim[1] - A_lim[0]) + A_lim[0]
        cube[2] = cube[2]*(Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]
        cube[3] = cube[3]*(offset_lim[1] - offset_lim[0]) + offset_lim[0]

    def log_likelihood(cube, ndim, nparams):
        L, d = np.random.choice(Ld_post, size=1)

        vals = models.offset_forward_model(r0, L, d, cube[0], w0, mu, cube[1], cube[2], coeff=cube[3])

        chisq = np.sum((vals - sig0)**2 / sig0_sd**2)

        return -chisq/2.0

    data = io.h5_2_dict(data_filename)
    Ld_post = np.vstack((Lpost, dpost)).T

    print(Ld_post.shape)

    ix0 = data['fit_ix']['0']

    r0 = data['r'][ix0]
    sig0 = data['sig'][ix0]
    sig0_sd = data['sig_sd'][ix0]

    F_lim = [17, 23]

    A_max = np.max(sig0)
    A_lim = [0.5*A_max, 2.0*A_max]

    Ti_lim = [0.025, 1.0]

    offset_lim = [0.0, 0.3]
    n_params = 4

    folder = path.abspath(output_folder)

    if test_plot:
        npts = 100
        nr = len(r0)

        vals = np.zeros((npts, nr))

        for i in range(npts):
            L, d = np.random.choice(Ld_post)

            cube = np.random.random(size=6)
            log_prior(cube, None, None)
            vals[i, :] = models.offset_forward_model(r0, L, d, cube[0], w0, mu, cube[1], cube[2], coeff=cube[3])

        fig, ax = plt.subplots()
        for i in range(npts):
            ax.plot(r0, vals[i, :], 'C0')
        ax.plot(r0, sig0, 'C1')
        plt.show()

    else:
        run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
            resume=resume, verbose=True, sampling_efficiency='q', n_live_points=75,
            outputfiles_basename=path.join(folder, 'full_'))
