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


def solve_Ar(shot_number, save_dir, Ld_directory, finesse_directory, order=0):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0] * (Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]
        cube[1] = cube[1] * (A_lim[1] - A_lim[0]) + A_lim[0]
        cube[2] = cube[2] * (V_lim[1] - V_lim[0]) + V_lim[0]
        #cube[3] = cube[3] * (Off_lim[1] - Off_lim[0]) + Off_lim[0] 

    def log_likelihood(cube, ndim, nparams):
        j = np.random.choice(i)
        L = Lpost[j]
        d = dpost[j]
        F = np.random.choice(Fpost)

        pred = model.forward3(rr, L, d, F, cube[0], mu, w, V=cube[2])*cube[1]
        #pred += cube[3]
        chisq = np.sum((ss-pred)**2 / ss_err**2)
        return -chisq / 2.0

    data_fname = "{0:07d}_data.json".format(shot_number)
    data_fname = join(save_dir, data_fname)
    with open(data_fname, 'r') as datafile:
        data = json.load(datafile, parse_float=np.float64)

    r = np.array(data['r'])
    s = np.array(data['sig'])
    print "arb subtraction"
    s -= 500.0
    idx = data['idx'][order]
    idx = range(idx[0]-50, idx[-1]+50) 

    rr = r[idx]
    ss = s[idx]
    ss_err = 0.03 * ss + 100.0

    A_lim = [0.25*np.max(ss), 3.0*np.max(ss)]
    Ti_lim = [0.025, 3.0]
    V_lim = [-10.0, 10.0]
    #Off_lim = [0.0, 1000.0]
    #print Off_lim
    Lpost, dpost = fp_helpers.read_Ld_results(Ld_directory)
    Fpost, _, _, _ = fp_helpers.read_finesse_results(finesse_directory)

    Lpost = Lpost[::30]
    dpost = dpost[::30]
    Fpost = Fpost[::10]
    i = range(len(Lpost))

    w = 487.98634
    mu = 39.948
    nparams = 3
    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=100,
            outputfiles_basename=join(save_dir,"fp_"), max_modes=500)



    #L = np.random.choice(Lpost)
    #d = np.random.choice(dpost)
    #F = np.random.choice(Fpost)
    #print L, d, F
    #pred = model.forward3(rr, L, d, F, 0.75, mu,  w, V=-0.15) * np.max(ss)*1.35
    ##with open(join(finesse_directory, 'calib_data.json'), 'r') as datafile:
    ##   calib_data = json.load(datafile, parse_float=np.float64)


    #fig, ax = plt.subplots() 
    #ax.plot(r, s)
    #ax.plot(rr, pred)
    ##ax1 = ax.twinx()
    ##ax1.plot(calib_data['r'], calib_data['sig'], 'g')
    #plt.show()
if __name__ == "__main__":
    solve_Ar(26248, "Calibration/Ti_saves/0026248_1/", "Calibration/Ld_saves/2017_09_25/",
            "Calibration/Finesse_saves/2017_09_25/", order=0)
