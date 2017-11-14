#!/usr/bin/env python
from __future__ import division
import pymultinest
import numpy as np
import os
import model
import time
import json
from os.path import join
import argparse
import random
import matplotlib.pyplot as plt
import fp_helpers

w = {'Th': 487.873302,
     'Ar': 487.98634}

mu = {'Th': 232.03806,
      'Ar': 39.948}

def solve_finesse_Ar(folder, data_fname, Ld_dir, order=0):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(F_lim[1]-F_lim[0]) + F_lim[0]
        cube[1] = cube[1]*(A_Ar_lim[1]-A_Ar_lim[0]) + A_Ar_lim[0]
        cube[2] = cube[2]*(A_Th_lim[1]-A_Th_lim[0]) + A_Th_lim[0]
        cube[3] = cube[3]*(Ar_Ti_lim[1]-Ar_Ti_lim[0]) + Ar_Ti_lim[0]

    def log_likelihood(cube, ndim, nparams):
        #L = np.random.choice(Lpost)
        #d = np.random.choice(dpost)
        j = np.random.choice(i)
        L = Lpost[j]
        d = dpost[j]

        pred = model.forward4(rr, L, d, cube[0], [Ti_Th, cube[3]], 
                [mu['Th'], mu['Ar']], [cube[1]*cube[2], cube[1]],
                [w['Th'], w['Ar']], V=0.0)
        chisq = np.sum((ss - pred)**2 / ss_err**2)
        #plt.plot(rr, ss, 'r')
        #plt.plot(rr, pred, 'b')
        #plt.show()
        return -chisq / 2.0


    with open(data_fname, 'r') as calibfile:
        calib_data = json.load(calibfile, parse_float=np.float64)
    r = np.array(calib_data['r'])
    sig = np.array(calib_data['sig'])
    
    idx = calib_data['idx'][order]
    rr = r[idx]
    ss = sig[idx]
    ss_err = 0.03*ss + 1.e3
    #ss_err = 0.03*ss + 1.e-1

    A_Ar_lim = [0.25*np.max(ss), 3.0*np.max(ss)]
    A_Th_lim = [0.01, 0.6]
    Ar_Ti_lim = [0.025, 2.0]
    F_lim = [10.0, 30.0]
    Ti_Th = 0.025 * 1000.0 / 300.0
    Lpost, dpost = fp_helpers.read_Ld_results(Ld_dir)
    Lpost = Lpost[::10]
    dpost = dpost[::10]
    i = range(len(Lpost))
    #log_likelihood([20.59, 0.1568e7, 0.2524, 0.2978], 0.0, 0.0)
    nparams = 4
    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=300,
            outputfiles_basename=join(folder,"fp_F_"), max_modes=500)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve for the etalon spacing d and camera focal length L.")
    parser.add_argument("finesse_dir", metavar="Finesse directory", type=str, 
            help="Name of finesse directory located in Calibration/Ld_saves/")
    parser.add_argument("Ld_dir", metavar="Ld directory", type=str, 
            help="Name of the Ld directory with MultiNest results.")
    parser.add_argument("--order", "-O", action='store', type=int, 
            help="Order number to solve for", default=0)
    args = parser.parse_args()


    finesse_dir = join("Calibration/Finesse_saves/", args.finesse_dir, "")
    Ld_dir = join("Calibration/Ld_saves/", args.Ld_dir, "")
    calib_fname = join(finesse_dir, "calib_data.json")

    solve_finesse_Ar(finesse_dir, calib_fname, Ld_dir, order=args.order)

