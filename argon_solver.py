#!/usr/bin/env python
import pymultinest
import numpy as np
import os
import multinest_plotting
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as mticker
#import seaborn.apionly as sns
import cPickle as pickle
import model
import plottingtools.core as ptools
import time
import json
import argparse
from calibration_solver import read_L_d_results, read_finesse_results, fix_save_directory
import random
from os.path import join

#def solve_Ar(savedir, data_fname, L, d, F):
#def solve_Ar(savedir, data_fname, etalon_dir, finesse_dir):
#def solve_Ar(savedir, data_fname, L, d, finesse_dir):
#def solve_Ar(savedir, data_fname, L, d):
def solve_Ar(savedir, data_fname, L, d, F):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]
        cube[1] = cube[1]*(V_lim[1] - V_lim[0]) + V_lim[0]
        cube[2] = cube[2]*(r0_lim[1] - r0_lim[0]) + r0_lim[0]
        cube[3] = cube[3]*(A_lim[1] - A_lim[0]) + A_lim[0]

    def log_likelihood(cube, ndim, nparams):

        #i = random.choice(n1)
        #j = random.choice(n2)

        #L = Ldpost[i, 0]
        #d = Ldpost[i, 1]
        #F = Fpost[j, 0]
        #F = 19.7714912064 
        linear_out = model.forward3(r, L, d, F, cube[0], muAr, wAr, nlambda=512, V=cube[1]) 
        linear_out *= cube[3]*np.exp(-(r / cube[2])**2)

        chisq = np.sum( (linear_out - s)**2 / s_sd**2 )
        return -chisq/2.0

    params = ["Ti", "V", "r0", "Amp"]
    labels = ["Ti (eV)", "V (km/s)", "r0 (px)", "Amp (Counts)"]
    prob_labels = ["P(Ti)", "P(V)", "P(r0)", "P(Amp)"]
    #shotnum = 9215

    param_info = {'params':params, 'labels': labels, 'prob labels': prob_labels}
    with open(savedir+"param_file.json", 'w') as param_config:
        json.dump(param_info, param_config, indent=4)

    with open(data_fname, 'r') as infile:
        data = json.load(infile, parse_float=np.float64)

    #with open(join(etalon_dir, "fp_Ld_post_equal_weights.dat"), 'r') as postfile:
    #    Ldpost = np.loadtxt(postfile, ndmin=2)

    #with open(join(finesse_dir, "fp_post_equal_weights.dat"), 'r') as postfile:
    #    Fpost = np.loadtxt(postfile, ndmin=2)
    
    #n1 = range(Ldpost.shape[0])
    #n2 = range(Fpost.shape[0])

    rarr = np.array(data["r"])
    sig = np.array(data["sig"])
    idx = np.array(data['idx'])
    r = rarr[idx]
    s = sig[idx]
    #s_sd = 0.01 * s + 100.0
    #s_sd = np.sqrt(s) + 10.0
    s_sd = 0.005 * s + 1.0
    #plt.plot(rarr, sig)
    #plt.plot(rarr[idx], sig[idx])
    #plt.show()
    # define limits for log_prior 
    #A_lim = [0.5 * data['A'], 5.*data['A']] # Counts
    A_lim = [0.01 * data['A'], 2*data['A']]
    Ti_lim = [0.025, 4.0]  # eV
    V_lim = [-10.0, 10.0]  # km/s
    r0_lim = [2000.0, 6000.0] # px
    nparams = len(params) 
    #print A_lim, sig.max()

    wAr = 487.98634
    muAr = 39.948
    #muAr = 40.0
    #plt.plot(rarr, sig)
    #plt.show()

    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=2000,
            outputfiles_basename=savedir+"fp_full_", max_modes=500)    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Argon Solver for the Fabry-Perot Interferometer.")
    parser.add_argument("--etalon-dir", action='store', type=str, default='saves/Ld_test7', dest='etalon_dir',
                        help="Directory of etalon solver from calibration_solver.py")
    parser.add_argument("--finesse-dir", action='store', type=str, default="saves/finesse_solver4",
                        dest="finesse_dir", help="Directory of finesse solver from calibration_solver.py")
    parser.add_argument("--savedir", "-s", action='store', type=str, default="solver_Ar1", dest='savedir',
                        help="Directory to store Ar solver files")
    args = parser.parse_args()

    #F = read_finesse_results(args.finesse_dir)
    savedir = fix_save_directory(args.savedir)
    etalon_dir = fix_save_directory(args.etalon_dir)
    finesse_dir = fix_save_directory(args.finesse_dir)
    #L, d = read_L_d_results(etalon_dir)
    #L = 148.599558491 / .004
    #d = 0.872203104696

    #L =148.715985377 / .004
    #d = 0.872201612198
    
    # old solution
    #L = 149.973834053 / .004
    #d = 0.877040367054
    #F = 22.2018789225
    ## new solution with bgfix
    #L = 149.975674078 / .004
    #d =0.877040359256
    #F = 21.9178512881

    L = 148.715306908/ .004
    d = 0.872201921125
    F = 20.0

    #L = 149.867 / .004
    #d = 0.8770405
    #F = 22.035
    #print "L (mm), d(mm), F", L*.004, d, F
    print "L (mm), d(mm)", L*.004, d
    #solve_Ar(savedir, "0009215_data.json", L, d, F)
    solve_Ar(savedir, "0009215_data2.json", L, d, F)
    #solve_Ar(savedir, "0009215_data.json", etalon_dir, finesse_dir)
    #solve_Ar(savedir, "0009215_data2.json", etalon_dir, finesse_dir)
    #solve_Ar(savedir, "0009215_data2.json", L, d, finesse_dir)
    #solve_Ar(savedir, "0009215_data2.json", L, d)
   # fname = join("synthetic_data/test1/", "Ar_noV_data.json")
    #fname = join("synthetic_data/test1/", "Ar_V_4.7_data.json")
    #fname = join("synthetic_data/test1/", "Ar_V_4.7_bgfix_data.json")
    #fname = join("synthetic_data/test1/", "Ar_V_4.7_Ti_1.1_data.json")
    #fname = join("synthetic_data/test1/", "Ar_V_4.7_Ti_1.1_smallnoise_data.json")
    #fname = join("synthetic_data/test1/", "Ar_V_n3.2_Ti_1.1_smallnoise_bgfix_data.json")
    #fname = join("synthetic_data/test1/", "Ar_noV_Ti_1.1_data.json")
    #fname = join("synthetic_data/test1/", "Ar_noV_Ti_1.1_smallnoise_bgfix_data.json")
    #fname = join("synthetic_data/test1/", "Ar_noV_Ti_1.1_smallnoise_data_bgfix.json")
    #solve_Ar(savedir, fname, L, d, F)

    # syn run 0: no V
    # syn run 1: V = 4.7 km/s 
    # syn run 2: same as 1 except with the correct F, also using the proper mu (no 40)
    # syn run 3: same as 1 using solved F, proper mu, changed error to 1% + 10.0
    # syn run 4: same V, but Ti=1.1 eV
    # syn run 5: no V, Ti=1.1eV 
    # syn run 6: back to V=4.7 and Ti=1.1 but with 2000 lives points intead of 500
    # syn run 7: same as 6 but with .05 noise scale instead of 2
    # syn run 8: ti=1.1 eV, V=-3.2 km/s
    # syn run 9: same as 8 with the correct finesse
    # syn run 10: same as 8 but with actual F, d, L
    # syn run 11: calib answers, run 8 data but with bgfix
    # syn run 12: 4.7 km/s , Ti=.58 (run 1), bgfix
    # syn run 13: same as 12, back to sqrt error
    # syn run 14:  same as 13, but with more wing in the fit 5% -> 1.5%
    # run 15: messing with error .005 * sig + 1.0
    # run 16: fixed warr in modeling
    # run 17: 4.7 km/s 0.58 eV with new calibration
    # run 18: -3.2 km/s 1.1 eV with new calibration
    # run 19: no V, 1.1 eV with new calibration

    #real data runs
    # run 25: 9215 data with combined L,d solution for real calib
