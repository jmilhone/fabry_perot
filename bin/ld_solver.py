from __future__ import division, absolute_import
import sys
sys.path.append("../")
import numpy as np
import argparse
from fabry.tools.file_io import h5_2_dict, dict_2_h5, prep_folder,read_Ld_results
from fabry.core.fitting import determine_fit_range, find_peak, find_maximum, gaussian_fit
from fabry.tools.plotting import ringsum_click, peak_plot, my_hist, tableau20_colors
from fabry.tools.images import get_data
import matplotlib.pyplot as plt
from fabry.core.models import peak_calculator
from os.path import join,isfile,abspath
from mpi4py import MPI
import pymultinest
import subprocess
from scipy.special import erf

px_size = 0.004# * 3

def get_pk_locations(r, s, s_sd, w, folder, pkerr=1, names=None, pkthresh=0.20, plotit=True):
    '''
    interactive clicking to get peak locations from ringsum

    Args:
        r (np.ndarray): r (pixel) array
        s (np.ndarray): ringsum signal array
        w (list of floats): list of wavelengths corresponding
                to peaks
        pkerr (float, default=0.2): error for peak locations in 
                units of percent of binsize at peak. Should be
                between 0.15 and 0.5.
        names (list of str, default=None): list of names of
                wavelengths, defaults to None which just 
                uses wavelength value as name
        pkthresh (float, default=0.1): threshold of maximum
                near guess to use for getting fit range for
                peak finding (see determine_fit_range help)
        plotit (bool, default=True): flag to plot the resulting
                peak fits
    Returns:
        peaks (dict): dictionary containg peak locations for
                all wavelengths in w, keys are the wavelengths
        peaks_sd (dict): dictionary containing peak standard deviations
                for all wavelengths in w, keys are the wavelengths
        orders (dict): dictionary containg ccd order numbers
                for all wavelengths in w, keys are the wavelengths
    '''
    binsizes = np.diff(np.concatenate(([0],r)))
    wnames = [str(x) for x in w]
    if names is None:
        names = wnames
    if type(pkthresh) not in [list, tuple]:
        pkthresh = [pkthresh]*len(w)

    peaks = {}
    orders = {}
    peaks_sd = {}
    for j,wstr in enumerate(wnames):
        pkguess,_ = ringsum_click(r**2, s,
            title='Click on the {0} peaks you wish to include'.format(names[j]))

        peak_list = []
        peaks_sd_list = []
        for i,pk in enumerate(pkguess):
            idxs = determine_fit_range(r,s,np.sqrt(pk),thres=pkthresh[i])
            #pk_r,pk_sd = find_peak(r[idxs],s[idxs],binsizes[idxs]/2.0,s_sd[idxs],plotit=True)
            #pk_r, pk_sd = gaussian_fit(r[idxs]**2, s[idxs])
            #print('kens', pk_r)
            #ix = np.abs(pk_r - r).argmin()
            pk_r = find_maximum(r[idxs], s[idxs], returnval=False)
            peak_list.append(pk_r)
            #peaks_sd_list.append(pk_sd)
            rvals = np.sort(r[idxs])
            ifind = np.searchsorted(rvals, pk_r, side='left')
            #print(rvals[ifind-1], pk_r, rvals[ifind])
            peaks_sd_list.append(2.0*(rvals[ifind] - rvals[ifind-1]))
            #print pk_r, pk_sd, binsizes[ix]
        peaks[wstr] = np.array(peak_list)
        peaks_sd[wstr] = np.array(peaks_sd_list)

        order_list = raw_input("Please enter the {0} peak orders you've just selected, separated by commas: ".format(names[j]))
        order_list = [int(x) for x in order_list.split(',')]
        orders[wstr] = order_list

    if plotit:
        name_dic = {}
        for i,ww in enumerate(wnames):
            name_dic[ww] = names[i]
        print(peaks)
        print(peaks_sd)
        print(orders)
        peak_plot(r, s, peaks, peaks_sd, orders)
    print(peaks_sd)
    return peaks, peaks_sd, orders

def ld_multinest_solver(peaks, peaks_sd, orders, basename, L_lim, d_lim, livepoints=1000, resume=True):

    ## one pixel is px_size mm so we need to convert the L_lim
    L_lim = [x/px_size for x in L_lim]
    ## d is computed in log space
    log_d_lim = [np.log10(x) for x in d_lim]

    wavelengths = peaks.keys()

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(L_lim[1]-L_lim[0]) + L_lim[0]
        #cube[1] = 10**(cube[1]*(log_d_lim[1] - log_d_lim[0]) + log_d_lim[0])
        cube[1] = cube[1]*(d_lim[1] - d_lim[0]) + d_lim[0]

    def log_likelihood(cube, ndim, nparams):
        chisq = 0.0
        prob = 1
        log_prob = 0.0
        for w in wavelengths:
            r = peak_calculator(cube[0], cube[1], float(w), orders[w])
            #chisq += np.sum( (r-peaks[w])**2 / (peaks_sd[w])**2 )
            for rj, sd, peak in zip(r, peaks_sd[w], peaks[w]):
        #        # prob += cumulative_distribution(rj-peak, 0.2*sd)
        #        #print(rj, peak, sd)
                #prob *= uniform_distribution(rj, peak, sd)
                prob *= tanh_distribution(rj, peak, sd)
                #log_prob += tanh_distribution(rj, peak, sd)
        #if prob <= 0.0:
        #    print('im here')
        #    return -np.inf

        #log_prob = np.log(prob)
        #return log_prob
        # print(prob)
        #print(prob)
        #if any(p <= 0 for p in prob):
        #    print 'im here'
        #    return -np.inf

        if prob <= 0.0:
            print('im here')
            return -np.inf
        else:
            return np.log(prob)

        #return sum(np.log(p) for p in prob)
        return -chisq / 2.0

    def cumulative_distribution(delta_x, sigma_x):
        return 0.5 * (1.0 + erf(-delta_x / (np.sqrt(2)*sigma_x)))

    def tanh_distribution(x, x0, sigma, plotit=False):
        left = x0 - sigma/2.0
        right = x0 + sigma/2.0
        weight = 0.1*sigma/2.0
        # print(left, right, x, sigma)
        prob = 0.5*(np.tanh((x-left)/weight) - np.tanh((x-right)/weight))
        if plotit:
            #print np.trapz(prob/sigma, x=x)
            #print(prob[0:10])
            #print(np.tanh(-10))
            fig, ax = plt.subplots()
            ax.plot(x, prob/sigma)
            plt.show()
        #return np.log((prob+1e-11) / sigma)
        #return np.log((prob / sigma) + 1e-15)
        return (prob+1e-9)/sigma

    def uniform_distribution(x, x0, sigma):
        if x0-sigma/2.0 < x < x0+sigma/2.0:
            return 1.0 / sigma
        else:
            return 1e-12 

    #k = '487.873302'
    #x0 = peaks[k][0]
    #sigma = peaks_sd[k][0]
    #sigma = 0.5
    #xarr = np.linspace(x0-3, x0+3, 1000)
    #print(peaks, peaks_sd)
    # _ = tanh_distribution(xarr, x0, sigma, plotit=True)

    nparams = 2
    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=resume, verbose=True, sampling_efficiency='model', n_live_points=livepoints,
            outputfiles_basename=basename, max_modes=500)

    # Lpost, dpost = read_Ld_results("/home/milhone/Research/python_FabryPerot/Data/2018_06_01/ArgonCalib/")
    # #idx = np.where(dpost > 0.8834)
    # #idx = np.where(Lpost > 154.65 / px_size)
    # #idx = np.where(Lpost < 154.45 / px_size)
    # idx = np.where(np.logical_and(Lpost > 154.45/px_size, Lpost < 154.65/px_size))
    # w = wavelengths[1]
    # #print(w)
    # order = 0
    # rvals = peak_calculator(Lpost, dpost, float(w), order)
    # peak_sd, peak = peaks_sd[w][order], peaks[w][order]
    # likeli = tanh_distribution(rvals, peak, peak_sd)
    # #print(np.min(np.exp(likeli)))
    # rarr = np.linspace(peak-2*peak_sd, peak+2*peak_sd, 1000)
    #fig, ax = plt.subplots()
    #ax.plot(rvals, np.exp(likeli), '.')
    #ax.plot(rvals[idx], np.exp(likeli)[idx], '.', color='C1')
    #ax.plot(rarr, np.exp(tanh_distribution(rarr, peak, peak_sd)), 'C2')
    #print(np.trapz(np.exp(tanh_distribution(rarr, peak, peak_sd)), rarr))
    #plt.show()
    #a=raw_input('press enter')

def ld_check(folder, bins=None, saveit=True):
    '''
    plots results from a Ld calibration

    Args:
        folder (str): folder that contains multinest data
        bins (int, default=20): number of bins to use in histogram plots
        saveit (bool, default=False): if true, will save
            plots in 'plots' subfolder
    '''
    data = h5_2_dict(join(folder, 'input_Ld_solver.h5'))
    Lpost, dpost = read_Ld_results(folder)
    if saveit:
        fig_folder = join(folder, 'Ld_solver_plots')
        prep_folder(fig_folder)

    hists = {}
    for w, pk in data['peaks'].items():
        h_list = []
        for i, n in enumerate(data['orders'][w]):
            h_list.append(peak_calculator(Lpost, dpost, float(w), n))
        hists[w] = h_list

    means = {}
    stds = {}
    for w, pk in data['peaks'].items():
        means_list = []
        stds_list = []
        for h in hists[w]:
            means_list.append(np.mean(h))
            stds_list.append(np.std(h))
        means[w] = means_list
        stds[w] = stds_list

    norder = max([len(x) for x in data['orders'].values()])
    nwaves = len(data['peaks'].keys())

    fig0, ax0 = plt.subplots(figsize=(10,6))
    peak_plot(data['r'],data['sig'],data['peaks'],data['peaks_sd'],data['orders'],fax=(fig0,ax0),anspks=means,anspks_sd=stds)

    plt.show(block=False)

    fig_, ax_ = plt.subplots()
    for i, w in enumerate(hists.keys()):
        for j, hist in enumerate(hists[w]):
            ax_.axvline(data['peaks'][w][j], zorder=15)
    ax_.plot(data['r'], data['sig'])
    plt.show(block=False)
    wtest = 468.335172
    # print wtest
    # print peak_calculator(0.382288362412689094E+05, 0.883875718242851827E+00, wtest, 0)
    # print peak_calculator(0.382288362412689094E+05, 0.883875718242851827E+00, wtest, 1)

    fig1, ax1 = plt.subplots(2,1,figsize=(6,8))
    my_hist(ax1[0], Lpost*px_size, bins=bins)
    ax1[0].set_xlabel('L (mm)', fontsize=18)
    ax1[0].set_ylabel('P (L)', fontsize=18)
    ax1[0].tick_params(labelsize=16)
    ax1[0].get_xaxis().get_major_formatter().set_useOffset(False)
    ax1[0].get_xaxis().get_major_formatter().set_scientific(False)
    my_hist(ax1[1], dpost, bins=bins)
    ax1[1].set_xlabel('d (mm)', fontsize=18)
    ax1[1].set_ylabel('P (d)', fontsize=18)
    ax1[1].tick_params(labelsize=16)
    ax1[1].get_xaxis().get_major_formatter().set_useOffset(False)
    ax1[1].get_xaxis().get_major_formatter().set_scientific(False)
    ax1[1].set_xticks(ax1[1].get_xticks()[::2])
    fig1.tight_layout()
    plt.show(block=False)

    fig2, ax2 = plt.subplots(norder,nwaves,figsize=(12,10))
    axx = ax2.reshape(norder,nwaves)
    for i, w in enumerate(hists.keys()):
        for j, hist in enumerate(hists[w]):
            my_hist(axx[j,i], hist, bins=bins)
            axx[j,i].axvline(data['peaks'][w][j], color='k',zorder=15)
            axx[j,i].axvspan(data['peaks'][w][j]-data['peaks_sd'][w][j]/2.0,data['peaks'][w][j]+data['peaks_sd'][w][j]/2.0,color='gray',alpha=0.4,zorder=15)
            axx[j,i].set_ylabel('Order {0:d}'.format(data['orders'][w][j]))
            axx[j,i].get_xaxis().get_major_formatter().set_useOffset(False)
            axx[j,i].get_xaxis().get_major_formatter().set_scientific(False)
            axx[j,i].tick_params(labelsize=16)
        axx[0,i].set_title('{0} nm'.format(w),fontsize=18)
        axx[-1,i].set_xlabel('R (px)', fontsize=18)
    fig2.tight_layout()

    #if saveit:
    #    fig0.savefig(join(fig_folder,'peaks.png'), dpi=400)
    #    fig1.savefig(join(fig_folder,'Ld_marginals.png'), dpi=400)
    #    fig2.savefig(join(fig_folder,'peak_histograms.png'), dpi=400)

    plt.show()

    #fig, ax = plt.subplots()
    ##hist = peak_calculator(Lpost, dpost, float(468.619458), 1)
    #hist = peak_calculator(Lpost, dpost, float(468.925186), 1)
    #my_hist(ax, hist)
    #plt.show()

if __name__ == "__main__":
    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    if rank == 0:
        parser = argparse.ArgumentParser(description='Performs L and d calibration.')
        parser.add_argument('folder', type=str, help='Folder where ringsum.h5 is from process_image. \
                This is where multinest outputs will be saved as well.')
        parser.add_argument('--wavelengths', '-w0', type=str, nargs='+', 
                default=['487.873302','487.98634'], help='wavelengths of peaks you\
                want to use to calibrate. Default is 487.873302 (ThI) and 487.98635 (ArII)')
        parser.add_argument('--pkerr', type=float, default=0.2,
                help='error for peak locations in units of percent of binsize at peak. \
                Should be between 1/6 and 1/2. Default is 1/5')
        parser.add_argument('--pkthresh', type=float, default=0.7,
                help='threshold for finding peak (can be a list corresponding to the various\
                peaks). Default is 0.4')
        parser.add_argument('--overwrite',action='store_true', help='allows you to overwrite\
                previously saved peak information')
        parser.add_argument('--L_lim', '-L', type=float, nargs=2, default=[145.,155.], help='limit\
                for fitting L in units of mm. Default is 145-155.')
        parser.add_argument('--d_lim', '-d', type=float, nargs=2, default=[0.87,0.89], help='limit\
                for fitting d in units of mm. Default is 0.87-0.89')
        parser.add_argument('--livepoints', type=int, default=2000, help='number of livepoints\
                for multinest to use. Default is 2000.')
        parser.add_argument('--no_solve',action='store_true', help='only writes peak information')
        parser.add_argument('--bins',type=int, default=20, help='number of bins to plot in histograms')
        args = parser.parse_args()

        fname = abspath(join(args.folder, 'input_Ld_solver.h5'))
        basename = abspath(join(args.folder, 'Ld_'))
        org_fname = abspath(join(args.folder, 'ringsum.h5'))

        if args.overwrite:
            resume = False
        else:
            resume = True
        resume = False 
        if not isfile(fname) or args.overwrite:
            if isfile(org_fname):
                data = h5_2_dict(org_fname)
                r = data['r']
                sig = data['sig']
                sig_sd = data['sig_sd']
                pk_dir = join(args.folder,'multinest_peaks/')
                prep_folder(pk_dir)
                peaks, peaks_sd, orders = get_pk_locations(r, sig, sig_sd, args.wavelengths, pk_dir, pkerr=args.pkerr, pkthresh=args.pkthresh)
                image_name = data['fname']
                dict_2_h5(fname, {'r':r, 'sig':sig, 'peaks':peaks, 'orders':orders, 'peaks_sd':peaks_sd})
            else:
                raise ValueError('{0} does not exist!'.format(org_fname))
        else:
           a = h5_2_dict(fname)
           peaks = a['peaks']
           orders = a['orders']
           peaks_sd = a['peaks_sd']

        inputs = {'L_lim':args.L_lim, 'd_lim':args.d_lim, 'livepoints':args.livepoints, 'pkerr':args.pkerr}
        dict_2_h5(fname, inputs, append=True)

        if args.no_solve:
            solver_in = None
        else:
            solver_in = {'peaks':peaks, 'orders':orders, 'basename':basename, 'peaks_sd':peaks_sd,
                    'L_lim':args.L_lim, 'd_lim':args.d_lim, 'livepoints':args.livepoints, 'resume':resume}
    else:
        solver_in = None

    solver_in = Comm.bcast(solver_in, root=0)
    if solver_in is not None:
        ld_multinest_solver(solver_in['peaks'], solver_in['peaks_sd'], solver_in['orders'],
                solver_in['basename'], solver_in['L_lim'], solver_in['d_lim'], 
                livepoints=solver_in['livepoints'], resume=solver_in['resume'])

    if rank == 0:
        ld_check(args.folder,bins=args.bins)
