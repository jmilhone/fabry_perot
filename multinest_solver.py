#!/usr/bin/env python
import pymultinest
import numpy as np
import os
import multinest_plotting
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker as mticker
import seaborn.apionly as sns
import cPickle as pickle
import model
import plottingtools.core as ptools
import time

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
#rcParams['xtick.major.pad'] = 12
err = 0.25
#          0    1     2           3       4         5        6        7
params = ['L', 'd', 'finesse', 'Ti_Th', 'Ti_Ar', 'Amp_Th', 'Amp_Ar', 'r0']
nparams = len(params)

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
    data_sd = ringsum[ind]
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
    
        rAr = [cube[0] * np.sqrt( m**2 / (np.floor(m) - 1.*j)**2 - 1.0) for j in range(npeaks)]
        rTh = [cube[0] * np.sqrt( n**2 / (np.floor(n) - 1.*j)**2 - 1.0) for j in range(npeaks)]
        chisq = sum( (rAr[j] - peaksAr[j])**2 / err**2 for j in range(npeaks) )
        chisq += sum( (rTh[j] - peaksTh[j])**2 / err**2 for j in range(npeaks) )

        return -chisq / 2.0
    savedir = fix_save_directory(savedir)

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

    L = marginals[0]['median']
    Lsd = marginals[0]['sigma']

    #d = None
    #dsd = None

    d = marginals[1]['median']
    dsd = marginals[1]['sigma']
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
    savedir = fix_save_directory("full_solver_run5")
    t00 = time.ctime(time.time())
    full_solver(savedir, "fp_ringsum_params.p")

    print "Time started: ", t00
    print "Time finished: ",time.ctime(time.time())
    # run 0 forgot to convert L to pixels
    # run 1 has the L fix, mixed up the wavelengths, mu's were right, should be fixed now
    # run 2
    # run 3 changed fall off from exp(-r/r0) to exp(-(r/r0)**2)  gaussian fall off
    # run 4 changed bin size from .2 pixels to .4
    # run 5 switched to model.forward2 where only one integral is done

