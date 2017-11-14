#!/usr/bin/env python
from __future__ import division
import pymultinest
import numpy as np
import os
import model
import time
import json
from check_caliibration_solution import plot_marginals
from os.path import join
import cPickle as pickle
import argparse
import random
import He_model2 as model2
import matplotlib.pyplot as plt
from scipy.stats import norm


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


def calculate_peaks(L, d, w, norders=4):
    m = 2.e6 * d / w
    m0 = np.floor(m)
    return [L * np.sqrt( m**2 / (m0 - j)**2 - 1.0) for j in range(norders)]
    #return [np.sqrt(2.0) * L * np.sqrt(1.0 - (m0-j) / m) for j in range(norders)]

#def solve_L_d(savedir, peak_file, d_lim=(0.88-0.01, 0.88+0.01),
def solve_L_d(savedir, peak_file, sigma_file, d_lim=(0.88-0.01, 0.88+0.01),
        L_lim=(145.0/.004, 155.0/.004)):

    def log_prior(cube, ndim, nparams):

        cube[0] = cube[0]*(L_lim[1]-L_lim[0]) + L_lim[0]
        cube[1] = 10**(cube[1]*(log_d_lim[1] - log_d_lim[0]) + log_d_lim[0])

    def log_likelihood(cube, ndim, nparams):
        chisq = 0.0
        for pk, w in zip(peak_labels, wavelengths):
            r = calculate_peaks(cube[0], cube[1], w, norders=npks)
            #chisq += sum((r[j] - peaks[pk][j])**2 / err**2 for j in range(npks))
            chisq += sum((r[j] - peaks[pk][j])**2 / sigma[pk][j]**2 for j in range(npks))

        return -chisq / 2.0

    log_d_lim = [np.log10(x) for x in d_lim]

    # changed to 2.0 for run 8, was at 1
    #err = 1.0 
    #err = 0.5

    params = ['L', 'd']
    labels = ["L (mm)", "d (mm)"]
    prob_labels = ["P(L)", "P(d)"]

    param_info = {'params':params, 'labels': labels, 'prob labels': prob_labels}
    with open(savedir+"param_file.json", 'w') as param_config:
        json.dump(param_info, param_config, indent=4)

    with open(peak_file, 'r') as infile:
        peaks = json.load(infile, parse_float=np.float64)
    with open(sigma_file, 'r') as infile:
        sigma = json.load(infile, parse_float=np.float64)
    npks = 3
    peak_labels = peaks.keys()
    wavelengths = [float(x) for x in peak_labels]
    n_params = 2
    pymultinest.run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=2000,
            outputfiles_basename=savedir+"fp_Ld_", max_modes=500)

def read_L_d_results(savedir, basename="fp_Ld_"):

    analyzer = pymultinest.Analyzer(n_params=2, outputfiles_basename=join(savedir,basename))

    stats = analyzer.get_mode_stats()
    local_log_ev = [x['local log-evidence'] for x in stats['modes']]
    ix = np.argmax(local_log_ev)

    mode = stats['modes'][ix]
    mode_vals = mode['maximum a posterior']

    L = mode_vals[0]
    d = mode_vals[1]
    print "L, d = ",L*.004, d
    return L, d


def read_finesse_results(savedir, basename="fp_", return_all=False):
    analyzer = pymultinest.Analyzer(n_params=2, outputfiles_basename=join(savedir,basename))

    stats = analyzer.get_mode_stats()
    local_log_ev = [x['local log-evidence'] for x in stats['modes']]
    ix = np.argmax(local_log_ev)

    mode = stats['modes'][ix]
    mode_vals = mode['maximum a posterior']

    F = mode_vals[0]
    if return_all:
        return mode_vals
    else:
        return F


def finesse_solver(savedir, param_fname, etalon_dir):#, L, d):
#def finesse_solver(savedir, param_fname, L, d):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(finesse_range[1] - finesse_range[0]) + finesse_range[0]
        cube[1] = cube[1]*(Ti_Ar_range[1] - Ti_Ar_range[0]) + Ti_Ar_range[0]
        cube[2] = cube[2]*(Th_amp_range[1] - Th_amp_range[0]) + Th_amp_range[0]
        cube[3] = cube[3]*(Ar_amp_range[1] - Ar_amp_range[0]) + Ar_amp_range[0]
        #cube[4] = cube[4]*(rscale_range[1] - rscale_range[0]) + rscale_range[0]

    def log_likelihood(cube, ndim, nparams):
        i = random.choice(nchoice)
        L = Ldpost[i, 0]
        d = Ldpost[i, 1]

        #rrr = np.zeros_like(rr)
        #for i in xrange(nr):
        #    rrr[i] = pdfs[i].rvs(1)[0]

        #linear_out = model.forward4(rrr, L, d, cube[0],
        #        [Ti_Th, cube[1]], [muTh, muAr], [cube[2], cube[3]], [wTh, wAr]) 
        #linear_out *= np.exp(-(rrr / cube[4])**2)

        linear_out = model.forward4(rr, L, d, cube[0],
                [Ti_Th, cube[1]], [muTh, muAr], [cube[2], cube[3]], [wTh, wAr]) 
        #linear_out *= np.exp(-(rr / cube[4])**2)
        chisq = np.sum( (linear_out - data)**2 / data_sd**2 )
        return -chisq/2.0 

    with open(join(etalon_dir, "fp_Ld_post_equal_weights.dat"), 'r') as postfile:
        Ldpost = np.loadtxt(postfile, ndmin=2)
    nchoice = range(Ldpost.shape[0])

    params = ['Finesse', 'Ti_Ar', 'Amp_Th', 'Amp_Ar', 'r0']
    labels = ["Finesse", r"$T_{i, Ar}$ (eV)", r"$A_{Th}$ (Counts)", 
        r"$A_{Ar}$ (Counts)", r"$r_0$ (px)"]
    prob_labels = ["P(Finesse)", r"P($T_{i,Ar}$)", r"P($A_{Th}$)", 
        r"P($A_{Ar}$)", r"P($r_0$)"]

    param_info = {'params':params, 'labels': labels, 'prob labels': prob_labels}
    with open(savedir+"param_file.json", 'w') as param_config:
        json.dump(param_info, param_config, indent=4)

    wTh = 487.873302
    wAr = 487.98634

    muAr = 39.948
    muTh = 232.03806

    with open(savedir+param_fname) as paramfile:
        param_dict = pickle.load(paramfile)

    finesse_range = param_dict["finesse_range"]
    Ti_Ar_range = param_dict["Ti_Ar_range"]
    Th_amp_range = param_dict["Th_amp_range"]
    Ar_amp_range = param_dict["Ar_amp_range"]
    rscale_range = param_dict["rscale_range"]
    rscale_range = [500.0, 3000.0]

    r = param_dict["binarr"]
    binsize = param_dict["binsize"]
    ringsum = param_dict["ringsum"]
    ind = param_dict["ind"]
    data = ringsum[ind]
    data_sd = 0.01 * data + 100.0
    rr = r[ind]
    bin_sigma = binsize[ind] / 6.0
    nr = len(rr)
    pdfs = []
    for i in xrange(nr):
        pdfs.append(norm(scale=bin_sigma[i], loc=rr[i]))

    Ti_Ar_range = [x*.025 / 300.0 for x in Ti_Ar_range]
    Ti_Th = 1000.0 * .025 / 300.0

    #nparams = len(params)
    nparams = len(params)-1
    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=100,
            outputfiles_basename=savedir+"fp_", max_modes=500)

def check_peaks(savedir, peak_file, basename="fp_Ld_"):
    with open(peak_file, 'r') as infile:
        peaks = json.load(infile, parse_float=np.float64)
    npks = 3
    peak_labels = peaks.keys()
    wavelengths = [float(x) for x in peak_labels]
    
    L, d = read_L_d_results(savedir)
    for pk, w in zip(peak_labels, wavelengths):
        print "\n{0} nm".format(pk)
        print peaks[pk]
        theo = calculate_peaks(L, d, w, norders=4)
        print theo
    print ""
    print calculate_peaks(L, d, 468.6194581)

def Th_two_image_solver(savedir, d_lim=(0.88-0.01, 0.88+0.01),
        L_lim=(147.0/.004, 153.0/.004)):

    def log_prior(cube, ndim, nparams):

        cube[0] = cube[0]*(L_lim[1]-L_lim[0]) + L_lim[0]
        cube[1] = 10**(cube[1]*(log_d_lim[1] - log_d_lim[0]) + log_d_lim[0])
        cube[2] = cube[2] * (Frange[1]-Frange[0]) + Frange[0]
        cube[3] = cube[3] * (A_Ar[1] - A_Ar[0]) + A_Ar[0]
        cube[4] = cube[4] * (A_Th[1] - A_Th[0]) + A_Th[0]
        cube[5] = cube[5] * (A_ThHe[1] - A_ThHe[0]) + A_ThHe[0]
        cube[6] = cube[6] * (4500.0 - 2000.0) + 2000.0
        cube[7] = cube[7] * (4500.0 - 2000.0) + 2000.0
        #cube[8] = cube[8] * (2.e6 - 0.0) + 0.0
        cube[8] = cube[8] * (1.25 - .025) + .025

    def log_likelihood(cube, ndim, nparams):
        chisq = 0.0
        QQ = (2.0 * cube[2] / np.pi)**2
        #Tii = .025 * 1000.0 / 300.0 
        #sig1 = 3.276569e-5 * np.sqrt(Tii/232.0) * 488.0
        #sig2 = 3.276569e-5 * np.sqrt(Tii/232.0) * 468.0 
        #norm1 = sig1 * np.sqrt(2.0 * np.pi)
        #norm2 = sig2 * np.sqrt(2.0 *  np.pi)
        #norm1 = 1.0
        #norm2 = 1.0
        model_Th = model.forward4(rr_Th, cube[0], cube[1],cube[2] ,
                #[Ti_Th, Ti_Ar], [muTh, muAr], [cube[4], cube[3]], [wTh, wAr]) 
                [Ti_Th, cube[8]], [muTh, muAr], [cube[4], cube[3]], [wTh, wAr]) 
        #model_Th += 6 * 4.1e6 / (1 + QQ) * norm1
        model_Th *= np.exp(-(rr_Th / cube[6])**2)

        model_ThHe = model.forward4(rr_ThHe, cube[0], cube[1], cube[2],
                [Ti_Th, ], [muTh, ], [cube[5], ], [468.6195, ])
        #model_ThHe += 9 * 3.6e6 / (1 + QQ) * norm2
        #model_ThHe += cube[8]
        model_ThHe *= np.exp(-(rr_ThHe / cube[7])**2)

        chisq += np.sum( (model_Th - data_Th)**2 / sd_Th**2)
        chisq += np.sum( (model_ThHe - data_ThHe)**2 / sd_ThHe**2)

        return -chisq / 2.0
    params = ['L', 'd', 'F', 'Amp_Ar', 'Amp_Th', 'Amp_ThHe', 'r0_Th' , 'r0_ThHe']
    labels = ['L', 'd', 'F', 'Amp_Ar', 'Amp_Th', 'Amp_ThHe', 'r0_Th' , 'r0_ThHe']
    #labels = ["Finesse", r"$T_{i, Ar}$ (eV)", r"$A_{Th}$ (Counts)", 
    #    r"$A_{Ar}$ (Counts)", r"$r_0$ (px)"]
    #prob_labels = ["P(se)", r"P($T_{i,Ar}$)", r"P($A_{Th}$)", 
    #    r"P($A_{Ar}$)", r"P($r_0$)"]
    prob_labels = ["P(L)", "P(d)", "P(F)", "P(A Ar)", "P (A Th)", "P (A ThHe)", "P(r0 Th)", "P(r0 ThHe)"]
    param_info = {'params':params, 'labels': labels, 'prob labels': prob_labels}

    with open(savedir+"param_file.json", 'w') as param_config:
        json.dump(param_info, param_config, indent=4)


    log_d_lim = [np.log10(x) for x in d_lim]
    Frange = [15.0, 30.0]

    wTh = 487.873302
    wAr = 487.98634

    muAr = 39.948
    muTh = 232.03806

    Ti_Th = 1000.0 * .025 / 300.0
    Ti_Ar = 0.20
    # synthetic
    #Ti_Ar = 1.0

    A_Ar = [0.5e7, 1e8]
    A_Th = [0.5e7, 1e8]
    A_ThHe = [1e6 ,2e7]
    # synthetic limits
    #A_Ar = [1e3, 20e3]
    #A_Th = [1e3, 20e3]
    #A_ThHe = [1e3, 20e3]

    #A_Ar = [0.01, 2.0]
    #A_Th = [0.01, 2.0]
    #A_ThHe = [0.01 ,2.0]
    #with open("Th_ThHe_data.p", 'rb') as infile:
    #with open(join("synthetic_data/test1/","Th_ThHe_data.p"), 'rb') as infile:
    #with open(join("synthetic_data/test1/","Th_ThHe_data_bgfix.p"), 'rb') as infile:
    #with open("Th_ThHe_data_bgfix_realdata.p", 'rb') as infile:
    #with open("Th_ThHe_data_bghack_realdata.p", 'rb') as infile:
    with open("Th_ThHe_data_bghack3_realdata.p", 'rb') as infile:
        data = pickle.load(infile)

    rr_ThHe = data['rThHe']
    data_ThHe = data['sThHe']
    #inds_ThHe = data['inds_ThHe']
    #plt.plot(rr_ThHe, data_ThHe)
    #plt.plot(rr_ThHe[inds_ThHe]**2, data_ThHe[inds_ThHe], 'c')
    #plt.show()
    #data_ThHe = data_ThHe[inds_ThHe]
    #rr_ThHe = rr_ThHe[inds_ThHe]

    #sd_ThHe = 0.01 * data_ThHe + 100.
    # run 10
    sd_ThHe = np.sqrt(1.e4 * data_ThHe) + 100.0
    # for synthetic data
    #sd_ThHe = np.sqrt(data_ThHe) + 10.0


    rr_Th = data['rTh']
    data_Th = data['sTh']
    #plt.plot(rr_Th, data_Th)
    #plt.show()
    #inds_Th = data['inds_Th']
    #data_Th = data_Th[inds_Th]
    #rr_Th = rr_Th[inds_Th]
    #sd_Th = 0.01 * data_Th + 100.
    # run 10
    sd_Th = np.sqrt(1.e4 * data_Th) + 100.0
    # for synthetic data
    #sd_Th = np.sqrt(data_Th) + 10.0
    #nparams = 8
    nparams = 9

    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=500,
            outputfiles_basename=savedir+"fp_", max_modes=500)




def He_Th_test_solver(savedir, d_lim=(0.88-0.01, 0.88+0.01),
        L_lim=(148.0/.004, 152.0/.004)):


    def log_prior(cube, ndim, nparams):

        cube[0] = cube[0]*(L_lim[1]-L_lim[0]) + L_lim[0]
        cube[1] = 10**(cube[1]*(log_d_lim[1] - log_d_lim[0]) + log_d_lim[0])
        cube[2] = cube[2] * (Frange[1]-Frange[0]) + Frange[0]
        cube[3] = cube[3] * (A_Ar[1] - A_Ar[0]) + A_Ar[0]
        cube[4] = cube[4] * (A_Th[1] - A_Th[0]) + A_Th[0]
        cube[5] = cube[5] * (A_He[1] - A_He[0]) + A_He[0]
    def log_likelihood(cube, ndim, nparams):
        chisq = 0.0

        model_Th = model.forward4(rr_Th, cube[0], cube[1],cube[2] ,
                [Ti_Th, Ti_Ar], [muTh, muAr], [cube[4], cube[3]], [wTh, wAr]) 

        chisq += np.sum( (model_Th - data_Th)**2 / sd_Th**2)

        model_He = model2.model_output(rr_He, cube[0], cube[1], Ti_He, He_amplitudes, Q=cube[2]) * cube[5]
        # Q is mislabeled

        chisq += np.sum( (model_He - data_He)**2 / sd_He**2)
        return -chisq / 2.0

    log_d_lim = [np.log10(x) for x in d_lim]

    wTh = 487.873302
    wAr = 487.98634

    muAr = 39.948
    muTh = 232.03806

    Ti_Th = 1000.0 * .025 / 300.0
    Ti_Ar = 0.20
    Ti_He = 0.857 #9174
    He_amplitudes = [0.5471, 0.12806, 1.4311, 2.2558] # 9174

    Frange = [15.0, 30.0]
    A_Ar = [1e7, 1e8]
    A_Th = [1e7, 1e8]
    A_He = [10000.0, 60000.0]

    with open("Th_He_data.p", "rb") as inputfile:
        data = pickle.load(inputfile)

    rr_He = data['rHe']
    data_He = data['sHe']
    inds_He = data['inds_He']
    data_He = data_He[inds_He]
    rr_He = rr_He[inds_He]
    sd_He = 0.05 * data_He + 100.

    rr_Th = data['rTh']
    data_Th = data['sTh']
    inds_Th = data['inds_Th']
    data_Th = data_Th[inds_Th]
    rr_Th = rr_Th[inds_Th]
    sd_Th = 0.01 * data_Th + 100.
    nparams = 6

    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=250,
            outputfiles_basename=savedir+"fp_", max_modes=500)



if __name__ == "__main__":
    #parser = argparse.ArgumentParser(description="Calibration Solver for the Fabry-Perot Interferometer")
    #group = parser.add_mutually_exclusive_group()
    #group.add_argument("--finesse", "-F", action='store_true',dest='finesse', 
    #        help="Run Finesse Solver.  Requires L and d solver to be done")
    #group.add_argument("--etalon", "-E", action='store_true', help="Run L and d solver.", dest='etalon')
    #parser.add_argument("--etalon-dir", action='store', type=str, 
    #        help="Save directory for etalon solver", default="Ld_test8", dest='etalon_dir')
    #parser.add_argument("--finesse-dir", action='store', type=str, dest='finesse_dir',
    #        help="Save directory for finesse solver", default="finesse_solver3")
    #parser.add_argument("--peak-fname", "-P", action='store', type=str, default='calibration_peaks5.json', 
    #        help="JSON file containing radial peak locations for each wavelength.", dest='peak_fname')
    #parser.add_argument("--param-fname", "-pm", action='store', dest='param_fname', type=str, 
    #        help="Params needed for finesse solver including data", default='fp_ringsum_params.p')
    #args = parser.parse_args()
    #savedir = fix_save_directory("ThAr_ThHe_testing_solver22")
    #Th_two_image_solver(savedir)
    #savedir = fix_save_directory("Ld_test20")
    #solve_L_d(savedir, "calibration_peaks8.json", "calibration_peaks8_sig.json")
    #savedir = fix_save_directory("finesse_solver10")
    #check_peaks(savedir, "calibration_peaks8.json", basename="fp_Ld_")
    #savedir = fix_save_directory("finesse_solver15")
    #finesse_solver(savedir, "fp_ringsum_params.p", "saves/Ld_test20/")
    #plot_marginals(savedir, basename="fp_Ld_", param_fname="param_file.json", save=False)
    #plot_marginals(savedir, basename="fp_", param_fname="param_file.json", save=False)
#    He_Th_test_solver(savedir)
    #if args.etalon:
    #    savedir = fix_save_directory(args.etalon_dir)
    #    solve_L_d(savedir, args.peak_fname)
    #    plot_marginals(savedir, basename="fp_Ld_", param_fname="param_file.json", save=False)
    #    check_peaks(savedir, args.peak_fname)
    #elif args.finesse:
    #    savedir = fix_save_directory(args.finesse_dir)
    #    etalon_dir = fix_save_directory(args.etalon_dir)
    #    L, d = read_L_d_results(etalon_dir)
    #    #finesse_solver(savedir, args.param_fname, etalon_dir)#, L, d)
    #    finesse_solver(savedir, args.param_fname, L, d)
    #else:
    #    print "You did not specify a calibration routine.  Exiting..."
    savedir = fix_save_directory("Ld_newtest0")
    peak_fname = "calibration_peaks_new0.json"
    peak_sig_fname = "calibration_peaks_new0_sig.json"
    solve_L_d(savedir, peak_fname, peak_sig_fname)
    check_peaks(savedir, peak_fname, basename="fp_Ld_")
    #solve_L_d(savedir, peak_fname)
    #L, d = read_L_d_results(savedir)
    #plot_marginals(savedir, basename="fp_Ld_", param_fname="param_file.json", save=False)
    
    #print "Predicted He peaks: ", calculate_peaks(L, d, 468.564736669)
    #print "Predicted Th,He peaks: ", calculate_peaks(L, d, 468.6195)
    #savedir = fix_save_directory("finesse_solver3")
    #finesse_solver(savedir, "fp_ringsum_params.p", L, d)

    # 7 has the two r0's
    # 8 has an offset for ThHe
    # 9 removed offset, messed with what portions are being fit in Ar
    # 10 going to use an error of sqrt(1.e4 * counts) + 100 instead of .01*counts + 100
    # 11 has 4 orders now
    # 12 using 500 lives points instead of 250
    # 13 synthetic data test1
    # 14 has wavelength fix and bg fix
    # 15 latest and greatest on old calibration
    # 16 offset Q correction terms
    # 17 trying to get Q correction to work properly
    # 18 bg hack
    # 19 added back Ti_Ar
    # 20 bg hack 2
    # 21 same as 20 but with sqrt error
    # 22 hack 3
