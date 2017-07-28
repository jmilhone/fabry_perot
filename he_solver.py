#!/usr/bin/env python
from __future__ import division
import pymultinest
import numpy as np
import He_model2 as model2
import json
import os
from os.path import join

def fix_save_directory(save_directory):
    save_dir = os.path.join("saves/", save_directory, "")
    if not os.path.isdir('saves'):
        try:
            os.mkdir('saves')
        except OSError, e:
            print "Error making directory, other process probably already maded it.  Moving on gracefully"

    if not os.path.isdir(save_dir):
        try:
            os.mkdir(save_dir)
        except OSError, e:
            print "Error making directory, other process probably already maded it.  Moving on gracefully"
    return save_dir



def solve_He_complex(savedir, datafname, L, d, F):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]
        cube[1] = cube[1]*(V_lim[1] - V_lim[0]) + V_lim[0]
        #cube[2] = cube[2]*(A_lim[1] - A_lim[0]) + A_lim[0]
        #cube[3] = cube[3]*(A1_lim[1] - A1_lim[0]) + A1_lim[0]
        #cube[4] = cube[4]*(A2_lim[1] - A2_lim[0]) + A2_lim[0]
        #cube[5] = cube[5]*(A3_lim[1] - A3_lim[0]) + A3_lim[0]
        #cube[6] = cube[6]*(A4_lim[1] - A4_lim[0]) + A4_lim[0]
        #cube[7] = cube[7]*(r0_lim[1] - r0_lim[0]) + r0_lim[0]
    
        cube[2] = 10**(cube[2]*(A_lim[1] - A_lim[0]) + A_lim[0])
        cube[3] = 10**(cube[3]*(A1_lim[1] - A1_lim[0]) + A1_lim[0])
        cube[4] = 10**(cube[4]*(A2_lim[1] - A2_lim[0]) + A2_lim[0])
        cube[5] = 10**(cube[5]*(A3_lim[1] - A3_lim[0]) + A3_lim[0])
        cube[6] = 10**(cube[6]*(A4_lim[1] - A4_lim[0]) + A4_lim[0])
    def log_likelihood(cube, ndim, nparams):
        #amplitudes = [cube[3], cube[4], cube[5], cube[6], 1.0]
        amplitudes = [cube[2], cube[3], cube[4], cube[5], cube[6]]
        #lin_out = cube[2]*model2.model_output(rr, L, d, cube[0], amplitudes, V=cube[1], F=F)
        lin_out = model2.model_output(rr, L, d, cube[0], amplitudes, V=cube[1], F=F)
        #lin_out *= np.exp(-(rr/cube[7])**2) 
        chisq = np.sum( (lin_out - ss)**2 / s_sd**2)
        return -chisq/2.0

    params = ["Ti", "V", "A0", "A1", "A2", "A3", "A4"]#,"r0"]
    labels = ["$T_i$ (eV)", "V (km/s)", "$A_0$ (Counts)","$A_1$ (Counts)","$A_2$ (Counts)",
            "$A_3$ (Counts)","$A_4$ (Counts)"]#,"$r_0$ (px)"]
    prob_labels = ["P($T_i$)", "P(V)","P($A_0$)","P($A_1$)","P($A_2$)","P($A_3$)",
            "P($4_0$)"]#,"P()",

    param_info = {'params':params, 'labels': labels, 'prob labels': prob_labels}
    with open(savedir+"param_file.json", 'w') as param_config:
        json.dump(param_info, param_config, indent=4)

    #Ti_lim = [0.025, 1.0]
    Ti_lim = [0.025, 2.0]
    V_lim = [-10.0, 10.0]
    #A_lim = [2., 5.]
    #A_lim = [2., 6.]
    A_lim = [-3.0, 3.0]
    A1_lim = [-3.0, 3.0]
    A2_lim = [-3.0, 3.0]
    A3_lim = [-3.0, 3.0]
    A4_lim = [-3.0, 3.0]

    #Ti_lim = [0.025, 4.0]
    #A_lim = [1000., 150000.]
    #A1_lim = [0.0, 10.0]
    #A2_lim = [0.0, 10.0]
    #A3_lim = [0.0, 10.0]
    #A4_lim = [0.0, 10.0]
    r0_lim = [1000.0, 5000.0]
    
    with open(datafname, 'r') as datafile:
        data = json.load(datafile, parse_float=np.float64)

    idx = data['idx']
    rr = np.array(data['r'])[idx]
    ss = np.array(data['sig'])[idx] 
    maxval = ss.max()
    s_sd = 0.02 * ss + 10.0 
    ss /= maxval
    s_sd /= maxval
    
    #nparams = 8
    nparams = 7
    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=2000,
            outputfiles_basename=savedir+"fp_")    

if __name__ == "__main__":

    L = 148.715306908/ .004
    d = 0.872201921125
    #F = 18.553482898
    # from old calibration Ar image only
    #F = 22.9272771929
    F = 20.0
    #datafname = "He_0009174_data.json"
    datafname = "He_0009174_bin1.0+_1order_data.json"
    #datafname = "He_0009174_bin1.0+_3order_data.json"
    savedir = fix_save_directory("He_solver_run9")
    solve_He_complex(savedir, datafname, L, d, F)


    # run 0: 9174 first attempt, 2% + 10.0
    # run 1: same as 0, but with F=22.93 from old Th calib solver (no He filter image)
    # run 2: evan's thesis has the last amplitude being the largest, going to normalize to that one instead of first 
    #   also fixed velocity error in intergral calculation
    #   didnt finish run 2
    # run 3: going to log space for rel amps, 10^{-3, +3} 
    #   didnt let finish
    # run 4: going to 500 lives points, also lowering lower limit for amplitudes, going to limit Ti to less than 1 eV
    # also going to make amplitude a log space as well
    # run 5: added a ind amp norm inside He_model2.py model_output function
    # run 6: changing finesse to 20.0
    # run 7: dividing data by max value, using 5 rel amps instead of 4 and 1.0
    # run 8, 1 order, 1.0 binsize instead of 2.0 for all prev runs  back to 2000 live points
    #   Ti can go back to 2.0 max, F=20.0
    # run 9, removing r0 because I'm doing one order
    # run 10, 3 ordres, r0 back in , 500 live points for speed
