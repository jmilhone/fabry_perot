#!/usr/bin/env python
from __future__ import division
import pymultinest
import numpy as np
import time
import json
from os.path import join
import argparse
import random
import matplotlib.pyplot as plt
from fp_helpers import peak_calculator

wavelengths = {'Th': 487.873302,
               'Ar': 487.98634,
               'He': 468.619458}

pk_error = 0.25


def Ld_solver(folder, peak_name):

    def log_prior(cube, ndim, nparams):

        cube[0] = cube[0]*(L_lim[1]-L_lim[0]) + L_lim[0]
        cube[1] = 10**(cube[1]*(log_d_lim[1] - log_d_lim[0]) + log_d_lim[0])

    def log_likelihood(cube, ndim, nparams):
        chisq = 0.0
        for w in wavelengths:
            r = peak_calculator(cube[0], cube[1], wavelengths[w], orders[w])
            chisq += np.sum( (r-peak_data[w])**2 / pk_error**2 )

        return -chisq / 2.0
    # BRB
    d_lim = (0.88-0.01, 0.88+0.01)
    L_lim = (145.0/0.004, 155.0/0.004)
    
    # PCX
    #d_lim = (0.4739-0.01, 0.4739+0.01)
    #L_lim = (100.0/0.004, 110.0/0.004)

    log_d_lim = [np.log10(x) for x in d_lim]


    with open(peak_name, 'r') as peakfile:
        peak_data = json.load(peakfile, parse_float=np.float64)

    Ar_peaks = peak_data['Ar']
    Th_peaks = peak_data['Th']
    He_peaks = peak_data['He']

    orders = {}
    orders['Ar'] = np.linspace(0.0, len(Ar_peaks)-1, len(Ar_peaks))
    orders['Th'] = np.linspace(0.0, len(Th_peaks)-1, len(Th_peaks))
    orders['He'] = np.linspace(0.0, len(He_peaks)-1, len(He_peaks))

    n_params = 2
    pymultinest.run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
            #resume=True, verbose=True, sampling_efficiency='model', n_live_points=2000,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=3000,
            outputfiles_basename=join(folder,"fp_Ld_"), max_modes=500)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve for the etalon spacing d and camera focal length L.")
    parser.add_argument("Ld_dir", metavar="Ld directory", type=str, 
            help="Name of Ld directory located in Calibration/Ld_saves/")
    args = parser.parse_args()

    pkname = "calibration_peaks.json"

    Ld_dir = join("Calibration/Ld_saves/", args.Ld_dir, "")
    pkname = join(Ld_dir, pkname)

    Ld_solver(Ld_dir, pkname)
