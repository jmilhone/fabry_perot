import pymultinest
import numpy as np
import os
import multinest_plotting
import matplotlib.pyplot as plt

err = 0.5

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

def L_d_results(analysis, peaksTh, peaksAr, wTh=487.873302, wAr=487.98634, plotit=True):
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
        names = ['L', 'd']
        p = multinest_plotting.PlotMarginalModes(analysis)
        post = analysis.get_equal_weighted_posterior()
        print post.shape
        #p.plot_marginal(0, with_ellipses=True, with_points=False, grid_points=100)
        #plt.show()
        for i in range(nparams):
            ax = plt.subplot(nparams, nparams, nparams*i + i +1)
            p.plot_marginal(i, with_ellipses=True, with_points=False, grid_points=75)
            ax.set_xlabel(names[i])
            ax.set_ylabel("Probability")
            #ax.plot(vals[i], 0, 'ro', ms=8)

            for j in range(i):
                ax = plt.subplot(nparams, nparams, nparams*j+i+1)
                cb = p.plot_conditional(i, j, with_ellipses=False, with_points=False, ax=ax)
                #ax.plot(vals[i], vals[j], 'ro', ms=8)
                ax.set_xlabel(names[i])
                ax.set_ylabel(names[j])

        plt.show()

    return L, d, Lsd, dsd
