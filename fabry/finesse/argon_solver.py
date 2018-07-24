from __future__ import print_function, division, absolute_import
import pymultinest
import numpy as np
from ..core.models import forward_model, offset_forward_model
import json
import argparse
import h5py
import matplotlib.pyplot as plt # For testing purposes only
from ..tools import file_io as io
import random
from os.path import abspath, join

w0 = (487.873302, 487.98634)
mu = (232.03806, 39.948)


def solver(output_folder, prior_filename, data_filename, Lpost, dpost, resume=True, test_plot=False):
    """
    MultiNest solver for point spread function calibration with the Argon filter

    Arguments:
        output_folder (str): path to folder containing input and output files for finesse solver
        prior_filename (str): path to file for the MultiNest prior (json format)
        data_filename (str): path to file for the input data 
        Lpost (np.ndarray): array of equally weighted marginal posterior L values
        dpost(np.ndarray): array of equally weighted marginal posterior d values
        resume (bool, optional): MultiNest will resume the calculation if True
        test_plot (bool, optional): for debugging purposes, allows the user to sample the prior and compare to input data
    """

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(F_lim[1] - F_lim[0]) + F_lim[0]
        cube[1] = cube[1]*(A_lim[1] - A_lim[0]) + A_lim[0]
        cube[2] = cube[2]*(Arel_lim[1] - Arel_lim[0]) + Arel_lim[0]
        cube[3] = cube[3]*(Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]

        for idx, (w, amp_lim) in enumerate(zip(w_extra, Arel_extra), 4):
            cube[idx] = cube[idx]*(amp_lim[1] - amp_lim[0]) + amp_lim[0]

    def log_likelihood(cube, ndim, nparams):
        # I want to fix this at some point
        # i = random.randint(0, nL-1)
        i = np.random.choice(nL)
        L = Lpost[i]
        d = dpost[i]
        # L = 0.380173301412519577E+05
        # d = 0.883628502371783142E+00
        # amps, w, mass, V, Ti = build_function_parameters(cube, nparams)

        amps = [cube[1]*cube[2], cube[1]]
        w = list(w0)
        mass = list(mu)
        Ti = [0.025, cube[3]]
        V = [0.0, 0.0]

        #vals = forward_model(r, L, d, cube[0], w, mass, amps, Ti,
        #        V, sm_ang=False, nlambda=2000)
        vals = offset_forward_model(r, L, d, cube[0], w, mass, amps, Ti,
                V, sm_ang=False, nlambda=2000)
        # trying to model offset here
        #vals += cube[1] * 0.15 / (1.0 + cube[0])

        chisq = np.sum((vals - sig)**2 / error**2)
        return -chisq / 2.0

    def build_function_parameters(cube, nparams):
        """
        Helper function for building some intermediate lists of parameters
        needed for the forward model.

        Note that you need to be careful with the cube parameter. It is not a
        python list! I believe it is some kind of fortran array. For example,
        you cannot call len() on it.
        """
        amps = [0.0 for _ in range(nparams-4+2)]
        amps[0] = cube[2]
        amps[1] = 1.0
        for idx, x in enumerate(amps[2:], 2):
            amps[idx] = cube[idx+2]
        #amps.extend([x for x in list(cube[4:])])
        amps = [x * cube[1] for x in amps]

        w = [x for x in w0]
        w += w_extra

        mass = [x for x in mu]
        mass += [mu[0] for _ in w_extra]

        V = [0.0 for _ in mass]

        #Ti = [0.025*1000.0/300.0, cube[3]]
        #Ti += [0.025*1000.0/300.0 for _ in w_extra]

        Ti = [0.025 for _ in w]
        Ti[1] = cube[3]

        return amps, w, mass, V, Ti

    with open(prior_filename, 'r') as infile:
        prior = json.load(infile, parse_float=np.float64)

    data = io.h5_2_dict(data_filename)

    nL = len(Lpost)
    ix = data['fit_ix'][0:-1:3]
    r = data['r'][ix]
    sig = data['sig'][ix]
    error = data['sig_sd'][ix]

    F_lim = prior['F_lim']
    A_lim = (0.6*np.max(sig), 1.4*np.max(sig))
    Arel_lim = prior['Arel_lim']
    Ti_lim = prior['Ti_lim']
    w_extra = prior['w_extra']
    Arel_extra = prior['Arel_extra']

    assert len(w_extra) == len(Arel_extra)

    n_params = 4 + len(w_extra)
    folder = abspath(output_folder)

    print('There are {0:d} paremeters for MultiNest'.format(n_params))

    if test_plot:
        npts = 30
        test_sig = np.zeros((npts, len(r)))
        for i in xrange(npts):
            j = random.randint(0, nL-1)
            L = Lpost[j]
            d = dpost[j]
            cube = [random.random() for _ in xrange(n_params)] 
            log_prior(cube, None, None)
            amps, w, mass, V, Ti = build_function_parameters(cube, n_params)
            test_sig[i, :] = forward_model(r, L, d, cube[0], w, mass, amps, Ti,
                                           V, sm_ang=False, nlambda=2000)

        fig, ax = plt.subplots()
        for i in xrange(npts):
            ax.plot(r, test_sig[i, :], 'C0')
        ax.errorbar(r, sig, yerr=error, fmt='', ecolor='C2', color='C1')
        plt.show()
    else:
        pymultinest.run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
                resume=resume, verbose=True, sampling_efficiency='model', n_live_points=100,
                outputfiles_basename=join(folder, 'finesse_'))


def full_solver(output_folder, prior_filename, data_filename, resume=True, test_plot=False):
    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(L_lim[1] - L_lim[0]) + L_lim[0]
        cube[1] = cube[1]*(d_lim[1] - d_lim[0]) + d_lim[0]
        cube[2] = cube[2]*(F_lim[1] - F_lim[0]) + F_lim[0]
        cube[3] = cube[3]*(A_lim[1] - A_lim[0]) + A_lim[0]
        cube[4] = cube[4]*(Arel_lim[1] - Arel_lim[0]) + Arel_lim[0]
        cube[5] = cube[5]*(Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]


    def log_likelihood(cube, ndim, nparams):
        #vals = forward_model(r, cube[0], cube[1], cube[2], w0, mu, [cube[3]*cube[4], cube[3]], [Ti_Th, cube[5]],
        #        [0.0, 0.0], sm_ang=False, nlambda=2000)
        vals = offset_forward_model(r, cube[0], cube[1], cube[2], w0, mu, [cube[3]*cube[4], cube[3]], [Ti_Th, cube[5]],
                [0.0, 0.0], sm_ang=False, nlambda=2000, coeff=0.5)

        chisq = np.sum((vals - sig)**2 / error**2)
        return -chisq / 2.0

    data = io.h5_2_dict(data_filename)

    ix = data['fit_ix'][0:-1:2]
    r = data['r'][ix]
    sig = data['sig'][ix]
    error = data['sig_sd'][ix]



    Ti_Th = 0.025

    L_lim = [135.0, 150.0]
    L_lim = [x / 0.004 for x in L_lim]

    d_lim = [0.86, 0.88]

    F_lim = [18.0, 22.0]

    Amax = np.max(sig)
    A_lim = [0.75*Amax, 2.0*Amax]

    Arel_lim = [0.3, 0.6]

    Ti_lim = [0.025, 0.3]

    n_params = 6
    folder = abspath(output_folder)

    if test_plot:
        npts = 100
        test_sig = np.zeros((npts, len(r)))
        for i in xrange(npts):
            cube = [random.random() for _ in xrange(n_params)] 
            log_prior(cube, None, None)
            test_sig[i, :] = forward_model(r, cube[0], cube[1], cube[2], w0, mu, [cube[3]*cube[4], cube[3]],
                                           [Ti_Th, cube[5]], [0.0, 0.0], sm_ang=False, nlambda=2000)

        fig, ax = plt.subplots()
        for i in xrange(npts):
            ax.plot(r, test_sig[i, :], 'C0')
        ax.errorbar(r, sig, yerr=error, fmt='', ecolor='C2', color='C1')
        plt.show()

    else:
        pymultinest.run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
                resume=resume, verbose=True, sampling_efficiency='model', n_live_points=75,
                outputfiles_basename=join(folder, 'full_'))








