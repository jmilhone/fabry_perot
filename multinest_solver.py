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
from check_caliibration_solution import plot_marginals
from scipy.stats import norm

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
#rcParams['xtick.major.pad'] = 12
err = 0.25
#          0    1     2           3       4         5        6        7
#params = ['L', 'd', 'finesse', 'Ti_Th', 'Ti_Ar', 'Amp_Th', 'Amp_Ar', 'r0']
#nparams = len(params)

Arparams = ['Ti', 'V', 'r0', 'Amp']
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

def solve_Ar(savedir, param_fname):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]
        cube[1] = cube[1]*(V_lim[1] - V_lim[0]) + V_lim[0]
        cube[2] = cube[2]*(r0_lim[1] - r0_lim[0]) + r0_lim[0]
        cube[3] = cube[3]*(A_lim[1] - A_lim[0]) + A_lim[0]

    def log_likelihood(cube, ndim, nparams):

        linear_out = model.forward3(r, L, d, F, cube[0], muAr, wAr, nlambda=512, V=cube[1]) 
        linear_out *= cube[3]*np.exp(-(r / cube[2])**2)

        chisq = np.sum( (linear_out - s)**2 / s_sd**2 )
        return -chisq/2.0

    params = ["Ti", "V", "r0", "Amp"]
    labels = ["Ti (eV)", "V (km/s)", "r0 (px)", "Amp (Counts)"]
    prob_labels = ["P(Ti)", "P(V)", "P(r0)", "P(Amp)"]
    shotnum = 9215
    param_info = {'params':params, 'labels': labels, 'prob labels': prob_labels}
    with open(savedir+"param_file.json", 'w') as param_config:
        json.dump(param_info, param_config, indent=4)

    #with open("0015676_data.json", 'r') as infile:
    with open("{0:07d}_data.json".format(shotnum), 'r') as infile:
        data = json.load(infile, parse_float=np.float64)

    L = data["L"]
    d = data["d"]
    F = data["F"]
    rarr = np.array(data["r"])
    sig = np.array(data["sig"])
    idx = np.array(data['idx'])
    r = rarr[idx]
    s = sig[idx]

    # run 0, 2
    #s_sd = .01 * s + 100.0  # 1% error, avoid division by zero
    # run 1, 3
    #s_sd = np.sqrt(np.abs(s))
    # run 4
    #s_sd = np.sqrt(3.5*np.abs(s))
    # run 5
    #s_sd = .03 * s + 100.0
    # run 6, 9
    s_sd = .01 * s + 100.0

    # define limits for log_prior 
    A_lim = [0.5 * data['A'], 5.*data['A']] # Counts
    Ti_lim = [0.025, 4.0]  # eV
    V_lim = [-10.0, 10.0]  # km/s
    r0_lim = [2000.0, 6000.0] # px
    nparams = len(params) 


    wAr = 487.98634
    muAr = 40.0


    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=2000,
            outputfiles_basename=savedir+"fp_full_", max_modes=500)    

def full_solver3(savedir, param_fname):


    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(L_range[1] - L_range[0]) + L_range[0]
        cube[1] = 10.0**(cube[1]*(np.log10(d_range[1]) - np.log10(d_range[0])) + np.log10(d_range[0]))
        cube[2] = cube[2]*(finesse_range[1] - finesse_range[0]) + finesse_range[0]
        cube[3] = cube[3]*(Ti_Ar_range[1] - Ti_Ar_range[0]) + Ti_Ar_range[0]
        cube[4] = cube[4]*(Th_amp_range[1] - Th_amp_range[0]) + Th_amp_range[0]
        cube[5] = cube[5]*(Ar_amp_range[1] - Ar_amp_range[0]) + Ar_amp_range[0]
        cube[6] = cube[6]*(rscale_range[1] - rscale_range[0]) + rscale_range[0]

    def log_likelihood(cube, ndim, nparams):
        # two integral approach
        #linear_out_Ar = model.forward(rr, cube[0], cube[1], cube[2], cube[4], muAr, wAr)
        #linear_out_Th = model.forward(rr, cube[0], cube[1], cube[2], cube[3], muTh, wTh)

        #linear_out_Ar *= cube[6] 
        #linear_out_Th *= cube[5]
        #linear_out = np.exp(-(rr / cube[7])**2 ) * (linear_out_Ar + linear_out_Th) 

        #chisq = np.sum( (linear_out - data)**2 / data_sd**2 )

        # single integral approach
        rrr = np.zeros_like(rr)
        for i in xrange(nr):
            rrr[i] = pdfs[i].rvs(1)[0]

        linear_out = model.forward4(rr, cube[0], cube[1], cube[2],
                [Ti_Th, cube[3]], [muTh, muAr], [cube[4], cube[5]], [wTh, wAr]) 
        linear_out *= np.exp(-(rr / cube[6])**2)
        chisq = np.sum( (linear_out - data)**2 / data_sd**2 )

        return -chisq/2.0

    params = ['L', 'd', 'Finesse', 'Ti_Ar', 'Amp_Th', 'Amp_Ar', 'r0']
    labels = ["L (mm)", "d (mm)", "Finesse", 
            r"$T_{i, Ar}$ (eV)", r"$A_{Th}$ (Counts)", r"$A_{Ar}$ (Counts)", r"$r_0$ (px)"]
    prob_labels = ["P(L)", "P(d)", "P(Finesse)", r"P($T_{i,Ar}$)", 
            r"P($A_{Th}$)", r"P($A_{Ar}$)", r"P($r_0$)"]

    param_info = {'params':params, 'labels': labels, 'prob labels': prob_labels}
    with open(savedir+"param_file.json", 'w') as param_config:
        json.dump(param_info, param_config, indent=4)

    wTh = 487.873302
    wAr = 487.98634
    
    #Changed for run3
    #muAr = 40.0
    #muTh = 232.0
    muAr = 39.948
    muTh = 232.03806

    with open(savedir+param_fname) as paramfile:
        param_dict = pickle.load(paramfile)

    finesse_range = param_dict["finesse_range"]
    Ti_Ar_range= param_dict["Ti_Ar_range"]
    Ti_Th= 1000.0 * .025 /300.0
    d_range= param_dict["d_range"]
    L_range = param_dict["L_range"]
    Th_amp_range= param_dict["Th_amp_range"]
    Ar_amp_range= param_dict["Ar_amp_range"]
    rscale_range= param_dict["rscale_range"]
    r = param_dict["binarr"]
    binsize = param_dict['binsize']
    ringsum = param_dict["ringsum"]
    ringsum_sd = param_dict["ringsum_sd"]
    ind = param_dict["ind"]
    # Convert from K to eV 
    Ti_Ar_range = [x*.025 / 300.0 for x in Ti_Ar_range]
    L_range = [x / 0.004 for x in L_range]

    data = ringsum[ind]
    #data_sd = ringsum_sd[ind]
    data_sd = 0.01 * data +100.0 #np.sqrt(3.5*np.abs(data))
    rr = r[ind]
    bin_sigma = binsize[ind] / 6.0
    nr = len(rr)
    pdfs = []
    for i in xrange(nr):
        pdfs.append(norm(scale=bin_sigma[i], loc=rr[i]))

    nparams = len(params)
    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=100,
            outputfiles_basename=savedir+"fp_full_", max_modes=500)


def full_solver2(savedir, param_fname):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(L_range[1] - L_range[0]) + L_range[0]
        cube[1] = 10.0**(cube[1]*(np.log10(d_range[1]) - np.log10(d_range[0])) + np.log10(d_range[0]))
        cube[2] = cube[2]*(finesse_range[1] - finesse_range[0]) + finesse_range[0]
        cube[3] = cube[3]*(Ti_Ar_range[1] - Ti_Ar_range[0]) + Ti_Ar_range[0]
        cube[4] = cube[4]*(Ar_amp_range[1] - Ar_amp_range[0]) + Ar_amp_range[0]


    def log_likelihood(cube, ndim, nparams):
        linear_out = model.forward2(rr, cube[0], cube[1], cube[2],
                [Ti_Th, cube[3]], [muTh, muAr], [cube[4]/rel_amp, cube[4]], [wTh, wAr]) 
        linear_out *= np.exp(-(rr / rscale)**2)
        chisq = np.sum( (linear_out - data)**2 / data_sd**2 )

        return -chisq/2.0

    params = ['L', 'd', 'Finesse', 'Ti_Ar', 'Amp_Ar',]
    labels = ["L (mm)", "d (mm)", "Finesse",  r"$T_{i, Ar}$ (eV)",  r"$A_{Ar}$ (Counts)"]
    prob_labels = ["P(L)", "P(d)", "P(Finesse)", r"P($T_{i,Ar}$)", r"P($A_{Ar}$)"]

    param_info = {'params':params, 'labels': labels, 'prob labels': prob_labels}
    with open(savedir+"param_file.json", 'w') as param_config:
        json.dump(param_info, param_config, indent=4)


    with open(savedir+param_fname) as paramfile:
        param_dict = pickle.load(paramfile)

    nparams = len(params)

    finesse_range = param_dict["finesse_range"]
    Ti_Ar_range= param_dict["Ti_Ar_range"] # in Kelvin
    Ti_Th= param_dict["Ti_Th"] # in eV, normally 1000 K * .025 eV / 300.0 K
    d_range= param_dict["d_range"]
    L_range = param_dict["L_range"]  # in mm 
    rel_amp = param_dict["rel_amp"]  # Th Amplitude = Ar amplitude / rel_amp
    Ar_amp_range= param_dict["Ar_amp_range"]
    rscale = param_dict["rscale"]
    r = param_dict["binarr"]
    ringsum = param_dict["ringsum"]
    ind = param_dict["ind"]
    
    rr = r[ind]
    data = ringsum[ind]
    data_sd = .03 * ringsum[ind] + 100.0


    Ti_Ar_range = [x*.025 / 300.0 for x in Ti_Ar_range] # to eV
    L_range = [x / 0.004 for x in L_range] # to px

    wTh = 487.873302
    wAr = 487.98634

    muAr = 40.0
    muTh = 232.0

    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=2000,
            outputfiles_basename=savedir+"fp_newfull_", max_modes=500)

def full_solver(savedir, param_fname):


    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(L_range[1] - L_range[0]) + L_range[0]
        cube[1] = 10.0**(cube[1]*(np.log10(d_range[1]) - np.log10(d_range[0])) + np.log10(d_range[0]))
        cube[2] = cube[2]*(finesse_range[1] - finesse_range[0]) + finesse_range[0]
        cube[3] = cube[3]*(Ti_Th_range[1] - Ti_Th_range[0]) + Ti_Th_range[0]
        cube[4] = cube[4]*(Ti_Ar_range[1] - Ti_Ar_range[0]) + Ti_Ar_range[0]
        cube[5] = cube[5]*(Th_amp_range[1] - Th_amp_range[0]) + Th_amp_range[0]
        cube[6] = cube[6]*(Ar_amp_range[1] - Ar_amp_range[0]) + Ar_amp_range[0]
        cube[7] = cube[7]*(rscale_range[1] - rscale_range[0]) + rscale_range[0]

    def log_likelihood(cube, ndim, nparams):
        # two integral approach
        #linear_out_Ar = model.forward(rr, cube[0], cube[1], cube[2], cube[4], muAr, wAr)
        #linear_out_Th = model.forward(rr, cube[0], cube[1], cube[2], cube[3], muTh, wTh)

        #linear_out_Ar *= cube[6] 
        #linear_out_Th *= cube[5]
        #linear_out = np.exp(-(rr / cube[7])**2 ) * (linear_out_Ar + linear_out_Th) 

        #chisq = np.sum( (linear_out - data)**2 / data_sd**2 )

        # single integral approach
        linear_out = model.forward2(rr, cube[0], cube[1], cube[2],
                [cube[3], cube[4]], [muTh, muAr], [cube[5], cube[6]], [wTh, wAr]) 
        linear_out *= np.exp(-(rr / cube[7])**2)
        chisq = np.sum( (linear_out - data)**2 / data_sd**2 )

        return -chisq/2.0

    params = ['L', 'd', 'Finesse', 'Ti_Th', 'Ti_Ar', 'Amp_Th', 'Amp_Ar', 'r0']
    labels = ["L (mm)", "d (mm)", "Finesse", r"$T_{i, Th}$ (K)", 
            r"$T_{i, Ar}$ (eV)", r"$A_{Th}$ (Counts)", r"$A_{Ar}$ (Counts)", r"$r_0$ (px)"]
    prob_labels = ["P(L)", "P(d)", "P(Finesse)", r"P($T_{i,Th}$)", r"P($T_{i,Ar}$)", 
            r"P($A_{Th}$)", r"P($A_{Ar}$)", r"P($r_0$)"]

    param_info = {'params':params, 'labels': labels, 'prob labels': prob_labels}
    with open(savedir+"param_file.json", 'w') as param_config:
        json.dump(param_info, param_config, indent=4)

    wTh = 487.873302
    wAr = 487.98634

    muAr = 40.0
    muTh = 232.0

    with open(savedir+param_fname) as paramfile:
        param_dict = pickle.load(paramfile)

    finesse_range = param_dict["finesse_range"]
    Ti_Ar_range= param_dict["Ti_Ar_range"]
    Ti_Th_range= param_dict["Ti_Th_range"]
    d_range= param_dict["d_range"]
    L_range = param_dict["L_range"]
    Th_amp_range= param_dict["Th_amp_range"]
    Ar_amp_range= param_dict["Ar_amp_range"]
    rscale_range= param_dict["rscale_range"]
    r = param_dict["binarr"]
    ringsum = param_dict["ringsum"]
    ringsum_sd = param_dict["ringsum_sd"]
    ind = param_dict["ind"]

    # Convert from K to eV 
    Ti_Ar_range = [x*.025 / 300.0 for x in Ti_Ar_range]
    Ti_Th_range = [x*.025 / 300.0 for x in Ti_Th_range]
    L_range = [x / 0.004 for x in L_range]

    data = ringsum[ind]
    #data_sd = ringsum_sd[ind]
    data_sd = 0.03 * data +100.0 #np.sqrt(3.5*np.abs(data))
    rr = r[ind]

    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=2000,
            outputfiles_basename=savedir+"fp_full_", max_modes=500)

def solve_L_d(savedir, peaksTh, peaksAr, d_lim=(0.88-0.01, 0.88+0.01),
        L_lim=(148.0/.004, 152.0/.004), wTh=487.873302, wAr=487.98634):

    def log_prior(cube, ndim, nparams):

        cube[0] = cube[0]*(L_lim[1]-L_lim[0]) + L_lim[0]
        cube[1] = 10**(cube[1]*(log_d_lim[1] - log_d_lim[0]) + log_d_lim[0])

    def log_likelihood(cube, ndim, nparams):


        #dm = [wAr / 2.e6 / cube[0] / ( 1.0 / np.sqrt(cube[0]**2 + peaksAr[j]**2) - 1.0 / np.sqrt(cube[0]**2 + peaksAr[j+1]**2)) for j in range(npeaks-1)]
        #dn = [wTh / 2.e6 / cube[0] / ( 1.0 / np.sqrt(cube[0]**2 + peaksTh[j]**2) - 1.0 / np.sqrt(cube[0]**2 + peaksTh[j+1]**2)) for j in range(npeaks-1)]
        #d = np.mean(dm + dn)
        #m = 2.e6 * d / wAr
        #n = 2.e6 * d / wTh

        m = 2.e6 * cube[1] / wAr
        n = 2.e6 * cube[1] / wTh
        k0 = 2.e6 * cube[1] / wUk0 
        k1 = 2.e6 * cube[1] / wUk1 

        rAr = [cube[0] * np.sqrt( m**2 / (np.floor(m) - 1.*j)**2 - 1.0) for j in range(npeaks)]
        rTh = [cube[0] * np.sqrt( n**2 / (np.floor(n) - 1.*j)**2 - 1.0) for j in range(npeaks)]
        rUk0 = [cube[0] * np.sqrt( k0**2 / (np.floor(k0) - 1.*j)**2 - 1.0) for j in range(npeaks)]
        rUk1 = [cube[0] * np.sqrt( k1**2 / (np.floor(k1) - 1.*j)**2 - 1.0) for j in range(npeaks)]

        chisq = sum( (rAr[j] - peaksAr[j])**2 / err**2 for j in range(npeaks) )
        chisq += sum( (rTh[j] - peaksTh[j])**2 / err**2 for j in range(npeaks) )
        chisq += sum( (rUk0[j] - peaksUk0[j])**2 / err**2 for j in range(npeaks) )
        chisq += sum( (rUk1[j] - peaksUk1[j])**2 / err**2 for j in range(npeaks) )

        return -chisq / 2.0

    params = ['L', 'd']
    labels = ["L (mm)", "d (mm)"]
    prob_labels = ["P(L)", "P(d)"]

    param_info = {'params':params, 'labels': labels, 'prob labels': prob_labels}
    with open(savedir+"param_file.json", 'w') as param_config:
        json.dump(param_info, param_config, indent=4)

    #savedir = fix_save_directory(savedir)
    wUk0 = 487.800942
    wUk1 = 467.3660927

    #peaksUk0 = np.sqrt(np.array([60022.8, 828886., 1.6008e6, 2.37548e6]))
    #peaksUk1 = np.sqrt(np.array([322803., 1.06251e6, 1.8048e6, 2.54998e6]))
    peaksUk0 = [244.52057396118195, 911.07036218754683, 1266.0274806441739, 1541.5990554793589] 
    peaksUk1 = [568.77586620923012, 1031.354410078168, 1343.9804619294503, 1597.2710700051123] 

    log_d_lim = [np.log10(x) for x in d_lim]
    npeaks = min(len(peaksTh), len(peaksAr))
    n_params = 2 

    print "Th peaks: ", peaksTh
    print "Ar peaks: ", peaksAr

    pymultinest.run(log_likelihood, log_prior, n_params, importance_nested_sampling=False,
            resume=True, verbose=True, sampling_efficiency='model', n_live_points=2000,
            outputfiles_basename=savedir+"fp_Ld_", max_modes=500)

    analysis = pymultinest.Analyzer(outputfiles_basename=savedir+"fp_Ld_", n_params=n_params)
    return analysis

def L_d_results(analysis, peaksTh, peaksAr, wTh=487.873302, wAr=487.98634, plotit=True, savedir=None):
    nparams = 2 
    stats = analysis.get_stats()
    marginals = stats['marginals']
    modes = stats['modes']
    print len(modes)
    local_log_ev = [x['local log-evidence'] for x in stats['modes']]
    ix = np.argmax(local_log_ev)

    #L = marginals[ix]['median']
    Lsd = marginals[0]['sigma']
    L = modes[ix]['maximum a posterior'][0]
    #d = None
    #dsd = None

    #d = marginals[1]['median']
    dsd = marginals[1]['sigma']
    d = modes[ix]['maximum a posterior'][1]
    #for x in marginals[1]:
    #    print x, marginals[1][x]
    print "L = {0} +/- {1}".format(L, Lsd) 
    print L * .004
    print "d = {0} +/- {1}".format(d, dsd) 
    npeaks = np.min((len(peaksTh), len(peaksAr)))
    #dm = np.array([wAr / 2.e6 / L / ( 1.0 / np.sqrt(L**2 + peaksAr[j]**2) - 1.0 / np.sqrt(L**2 + peaksAr[j+1]**2)) for j in range(npeaks-1)])
    #dn = np.array([wTh / 2.e6 / L / ( 1.0 / np.sqrt(L**2 + peaksTh[j]**2) - 1.0 / np.sqrt(L**2 + peaksTh[j+1]**2)) for j in range(npeaks-1)])
    #d = np.mean(dm + dn)
    #m = 2.e6 * dm / wAr
    #n = 2.e6 * dn / wTh
    m = 2.e6 * d / wAr
    n = 2.e6 * d / wTh
    print " "
    print "Ar d: ", d
    print "Th d: ", d
    print " "
    print "Ar m: ", m
    print "Th n: ", n
    print " " 


    #d = 0.5 * (np.mean(dm) + np.mean(dn))  # dm and dn are the same length
    #d = (np.sum(dm) + np.sum(dn)) / (len(dm) + len(dn))
    #print d
    mAr =  2.e6 * d / wAr
    nTh = 2.e6 * d / wTh
    

    Arpks = [L * np.sqrt(mAr**2 / (np.floor(mAr)-j)**2 - 1.0) for j in range(npeaks)] 
    Thpks = [L * np.sqrt(nTh**2 / (np.floor(nTh)-j)**2 - 1.0) for j in range(npeaks)] 
    print "Ar comparison"
    print "Sol: ",Arpks
    print "Act: ",peaksAr
    print "m: ",mAr
    print ""
    print "Th comparison"
    print "Sol: ",Thpks
    print "Act: ",peaksTh
    print "m: ",nTh

    chisqAr = np.sum( (Arpks[j] - peaksAr[j])**2 / err**2 for j in range(npeaks) )
    chisqTh = np.sum( (Thpks[j] - peaksTh[j])**2 / err**2 for j in range(npeaks) )
    chisq = chisqAr + chisqTh

    print chisqAr, chisqTh, chisq
    
    if plotit:
        names = ['L (px)', 'd (mm)']
        p = multinest_plotting.PlotMarginalModes(analysis)
        post = analysis.get_equal_weighted_posterior()
        Lpost = post[:, 0]
        dpost = post[:, 1]
        plt.plot(dpost)
        plt.show()

        print post.shape
        #p.plot_marginal(0, with_ellipses=True, with_points=False, grid_points=100)
        #plt.show()
        for i in range(nparams):
            ax = plt.subplot(nparams, nparams, nparams*i + i +1)
            #p.plot_marginal(i, with_ellipses=True, with_points=False, grid_points=75)
            if i==0:
                a, b, _= ax.hist(Lpost, bins=40, normed=True)
                #sns.distplot(Lpost, bins=40, norm_hist=True, kde=False)

                #LL = np.linspace(L-3*Lsd, L+3.*Lsd, 100)
                #P_L = np.exp(-0.5 * (LL-L)**2 / Lsd**2 ) / Lsd / np.sqrt(2.0 * np.pi)
                #ax.plot(LL, P_L, '--r')
            else:
                a, b, _ = ax.hist(dpost, bins=40, normed=True)
                #temp = sns.distplot(dpost, bins=40, norm_hist=True, kde=False)
                #newfig, newax = plt.subplots()
                #newax.plot(x,y)
                
                #dd = np.linspace(d-3.*dsd, d+3.*dsd, 100)
                #P_d = np.exp(-0.5 * (dd-d)**2 / dsd**2 ) / dsd / np.sqrt(2.0 * np.pi)
                #ax.plot(dd, P_d, '--r')
            ax.set_xlabel(names[i], fontsize=16)
            ax.set_ylabel("Probability", fontsize=16)
            ax.tick_params(labelsize=16)
            #ax.plot(vals[i], 0, 'ro', ms=8)
            #if i==1:
            #    ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
            if i==0:
                ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
            myLocator = mticker.MaxNLocator(3)
            ax.xaxis.set_major_locator(myLocator)

            for j in range(i):
                ax = plt.subplot(nparams, nparams, nparams*j+i+1)
                #cb = p.plot_conditional(i, j, with_ellipses=False, with_points=False, ax=ax)
                _, _, _, im = plt.hist2d(dpost, Lpost, bins=75, normed=True)
                #g = sns.JointGrid(dpost, Lpost)
                #g.plot_joint(sns.kdeplot)
                #print type(im)
                cb = plt.colorbar(im )
                cb.set_label("Probability", fontsize=16)
                #ax.plot(vals[i], vals[j], 'ro', ms=8)
                ax.set_xlabel(names[i], fontsize=16)
                ax.set_ylabel(names[j], fontsize=16)
                #ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
                myLocator = mticker.MaxNLocator(3)
                ax.xaxis.set_major_locator(myLocator)
                ax.tick_params(labelsize=16)
        plt.tight_layout()
        if savedir:
            print "Going to save figure"
            save_dir = fix_save_directory(savedir)
            fig = plt.gcf()
            fig.set_size_inches(8,8)
            fig.savefig(save_dir + "marginals.pdf")
        plt.show()
    fig, ax = plt.subplots()
    ax.hist(dpost, bins=100, normed=True)
    ax.set_xlabel("d (mm)", fontsize=18)
    ax.set_ylabel("P(d)", fontsize=18)
    ptools.add_thick_box(ax, minor=False) 
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.tight_layout()
    fig.savefig("d_100bins.pdf")
    plt.show()
    return L, d, Lsd, dsd

if __name__ == "__main__":
    #savedir = fix_save_directory("full_solver_run18")
    t00 = time.ctime(time.time())
    #peaksTh = np.array([644.383, 1089.1, 1399.2, 1653. ]) 
    #peaksAr = np.array([733.8, 1144.7, 1443.5, 1691. ]) 
    #peaksTh = [644.45936512469882, 1089.370588979835, 1399.4961244393533, 1653.9688234764596]
    #peaksAr = [734.15309407063455, 1145.004498782735, 1443.8927124858421, 1691.3157088640837] 
    #savedir = fix_save_directory("Ld_test4")
    #analysis = solve_L_d(savedir, peaksTh, peaksAr)
    #L_d_results(analysis, peaksTh, peaksAr)
    #plot_marginals(savedir, basename="fp_Ld_", param_fname="param_file.json", save=False)
    #savedir = fix_save_directory("solveAr_0")
    #solve_Ar(savedir, None)
    #savedir = fix_save_directory("full_solver_run17")
    #full_solver(savedir, "fp_ringsum_params.p")

    #savedir = fix_save_directory("new_full_solver_run0")
    #full_solver2(savedir, "fp_ringsum_params.p")
    #savedir = fix_save_directory("solver3_run16")
    savedir = fix_save_directory("NewCalibration_Ar_run0")
    full_solver3(savedir, "fp_ringsum_params.p")
    #savedir = fix_save_directory("new_full_solver_run5")
    #full_solver2(savedir, "fp_ringsum_params.p")
    print "Time started: ", t00
    print "Time finished: ",time.ctime(time.time())
    # run 0 forgot to convert L to pixels
    # run 1 has the L fix, mixed up the wavelengths, mu's were right, should be fixed now
    # run 2
    # run 3 changed fall off from exp(-r/r0) to exp(-(r/r0)**2)  gaussian fall off
    # run 4 changed bin size from .2 pixels to .4
    # run 5 switched to model.forward2 where only one integral is done
    # run 6 changed bin size to 1.0 pixels
    # run 7 subtracted min value from ring sum
    # run 8 increased posterior of L to 152 mm, only did the first 3 orders
    # run 9, 3 orders, 1.0 pixels, no offset subtraction
    # run 10, synthetic data
    # run 11, same synthetic keeping all data points inbeween argon and th peaks
    # run 12, trying run 11 on real thorium lamp image
    # run 13, 12 was really messed up.  I changed the prior for the amplitudes from .5 to 5 to .5 to 2 times the amplituded guess
    # run 14, running 13 with an offset subtracted
    # run 15, same as run 14 except limiting F from 10 to 30, and r0 up to 5000
    # I found an type in the sd calculation
    # run 16 running with the suggested noise factor 
    # run17 was having trouble getting 16 to converge, trying 1% this time (unlike earlier runs)
    # run18, 17 is still going.  going to run it with 3% error

    # Ar solving
    # run 0: used 1% error with an offset of 100 counts
    # run 1: going to use the sqrt(counts) as the error
    # run 2: going back to 1% error, going to increase the prior for the amplitude
    # run 3: going back to sqrt(counts) with the increase in amplitude prior in run 2
    # run 4: going to use a noise factor of 3.5 aka error = sqrt(3.5 * counts)
    # run 5: going to use 1% error with full solver run17 solution
    # run 6: using 3% error to redo run 5
    # run 7: made a modification to the limits for the Amplitude
    # run 8: changed bottom finder from 2.5% to 10% in hopes to get the peaks fit better
    # run 9: running analysis on shotnum 9214
    # run 10: running analysis on shotnum 9215
    # run 11: running analysis on shotnum 9756
    # run 12: using proper ringsum r bins for shotnum 9756
    # run 13: changed the minimum subtraction to the last point in the array
    # run 14: running against forward model with forward model calibration solve as well
    # run 15: same Ti but with V=3.4 km/s
    # run 16: trying it with a negative velocity and small Ti
    # run 17: same as 16 but with perfect L and d input
    # run 18: shot=4, didnt subtract offset, using full solver solution
    # run 19: shot=5, went back to constant area ring sum
    # run 20: shot=6, no offset subraction, using solver3_run4 with no offset subtraction
    # run 21: shot=7, using solver3_run5
        # Time started:  Thu May 25 09:39:36 2017
        # Time finished:  Thu May 25 09:44:26 2017
    # run 22, shot=7, using solver3_run5, checking speed against foward3 vs forward in model.py
        # Time started:  Thu May 25 11:33:08 2017
        # Time finished:  Thu May 25 11:37:27 2017
        # 22 was 31 seconds faster
    # run 23, shot=8, looking at using 5% for thres2 instead of 10%

    # 5 param solver
    # run 0: initial run with 3% error


    # 7 param solver (no Ti_Th)
    # run0: initial run with 1% error
        # 16 cores on dave
        #Time started:  Mon May 15 08:26:59 2017
        #Time finished:  Mon May 15 22:22:40 2017
    # runn1: running on old thorium calibration data
        # 22 cores on dave
    # run2: trying it against the forward model
    # run3: reducing the prior size
    # run4: rerunning 3 but without subtracting an offset
    # run14: running on calibration from 9114
    # run15: Trying the small angle approx in forward4
    # run16: Added a 4th order

