from __future__ import print_function, division, absolute_import
import pymultinest
import numpy as np
from ..core.models import forward_model, peak_calculator
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

    Arguments:
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
    sig = data['sig'][ix]
    error = data['sig_sd'][ix]

    w_extra = prior['w_extra']

    n_params = 4 + len(w_extra)

    analyzer = pymultinest.Analyzer(n_params, outputfiles_basename=join(finesse_folder, "finesse_"))
    mode_stats = analyzer.get_mode_stats()

    # Find the mode with the lowest local log-evidence
    minval = np.inf
    ix = 0
    for idx, stat in enumerate(mode_stats['modes']):
        if stat['local log-evidence'] < minval:
            print('im here')
            ix = idx
            minval = stat['local log-evidence']

    mode = mode_stats['modes'][ix]['mean']

    amps = np.zeros(2 + len(w_extra))
    amps[0] = mode[2]
    amps[1] = 1.0
    amps[2:] = mode[4:]
    amps *= mode[1]
    amps = list(amps)

    w = [x for x in w0]
    w += w_extra

    mass = [x for x in mu]
    mass += [mu[0] for _ in w_extra]

    V = [0.0 for _ in mass]

    Ti = [0.025*1000.0/300.0, mode[3]]
    Ti += [0.025*1000.0/300.0 for _ in w_extra]

    L = 0.380173301412519577E+05
    d = 0.883628502371783142E+00

    vals = forward_model(r, L, d, mode[0], w, mass, amps, Ti, V, sm_ang=False, nlambda=2000)
    
    
    fig, ax = plt.subplots()
    ax.errorbar(r, sig, yerr=error, ecolor='C2', color='C0', label='Data')
    ax.plot(r, vals, color='C1', label='Fit')
    styles = ['-', '--', '-.']
    for idx, wavelength in enumerate(w):
        n, m = divmod(idx, 10)
        if wavelength in [488.185311, 488.32348]:
            order = 1
        else:
            order = 0
        ax.axvline(peak_calculator(L, d, wavelength, order), label=str(wavelength), color='C{0:d}'.format(m), linestyle=styles[n])
    ax.legend()
    plt.show(block=False)

    fig, ax = plt.subplots()
    ax.plot(r, (vals - sig)**2 / error**2)
    plt.show()

