import numpy as np
import argparse
from tools.file_io import h5_2_dict, dict_2_h5, prep_folder,read_Ld_results
from core.fitting import determine_fit_range, find_maximum
from tools.plotting import ringsum_click, peak_plot, my_hist, tableau20_colors
from tools.images import get_data
import matplotlib.pyplot as plt
from core.models import peak_calculator
from os.path import join,isfile,abspath
from mpi4py import MPI
import pymultinest
import subprocess

def get_pk_locations(r, s, w, names=None, pkthresh=0.10, plotit=True):
    '''
    interactive clicking to get peak locations from ringsum

    Args:
        r (np.ndarray): r (pixel) array
        s (np.ndarray): ringsum signal array
        w (list of floats): list of wavelengths corresponding
                to peaks
        names (list of str, default=None): list of names of
                wavelengths, defaults to None which just 
                uses wavelength value as name
    Returns:
        peaks (dict): dictionary containg peak locations for
                all wavelengths in w, keys are the wavelengths
        orders (dict): dictionary containg ccd order numbers
                for all wavelengths in w, keys are the wavelengths
    '''
    wnames = [str(x) for x in w]
    if names is None:
        names = wnames
    if type(pkthresh) not in [list, tuple]:
        pkthresh = [pkthresh]*len(w)

    peaks = {}
    orders = {}
    for i,wstr in enumerate(wnames):
        pkguess,_ = ringsum_click(r**2, s, 
                title='Click on {0} peaks you wish to include'.format(names[i]))
        
        peak_list = []
        for pk in pkguess:
            idxs = determine_fit_range(r,s,np.sqrt(pk),thres=pkthresh[i])
            pk_r = find_maximum(r[idxs],s[idxs])
            peak_list.append(pk_r)
        peaks[wstr] = np.array(peak_list)
        
        order_list = raw_input("Please enter the {0} peak orders you've just selected, separated by commas: ".format(names[i]))
        order_list = [int(x) for x in order_list.split(',')]
        orders[wstr] = order_list

    if plotit:
        name_dic = {}
        for i,ww in enumerate(wnames):
            name_dic[ww] = names[i]
        peak_plot(r, s, peaks, orders)
    
    return peaks, orders

def ld_multinest_solver(peaks, orders, basename, pk_error, L_lim, d_lim, livepoints=1000, resume=True):
    
    ## one pixel is 0.004 mm so we need to convert the L_lim
    L_lim = [x/0.004 for x in L_lim]
    ## d is computed in log space
    log_d_lim = [np.log10(x) for x in d_lim]

    wavelengths = peaks.keys()

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(L_lim[1]-L_lim[0]) + L_lim[0]
        cube[1] = 10**(cube[1]*(log_d_lim[1] - log_d_lim[0]) + log_d_lim[0])

    def log_likelihood(cube, ndim, nparams):
        chisq = 0.0
        for w in wavelengths:
            r = peak_calculator(cube[0], cube[1], float(w), orders[w])
            chisq += np.sum( (r-peaks[w])**2 / pk_error[w]**2 )
        return -chisq / 2.0
    
    nparams = 2
    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=resume, verbose=True, sampling_efficiency='model', n_live_points=livepoints,
            outputfiles_basename=basename, max_modes=500)

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

    norder = max([len(x) for x in data['orders'].values()])
    nwaves = len(data['peaks'].keys())

    fig0, ax0 = plt.subplots(figsize=(10,6))
    peak_plot(data['r'],data['sig'],data['peaks'],data['orders'],fax=(fig0,ax0),pkerr=data['pkerr'])
    plt.show(block=False)

    fig1, ax1 = plt.subplots(2,1,figsize=(6,8))
    my_hist(ax1[0], Lpost*0.004, bins=bins)
    ax1[0].set_xlabel('L (mm)')
    ax1[0].set_ylabel('P (L)')
    ax1[0].get_xaxis().get_major_formatter().set_useOffset(False)
    ax1[0].get_xaxis().get_major_formatter().set_scientific(False)
    ax1[0].set_xlim(data['L_lim'])
    my_hist(ax1[1], dpost, bins=bins)
    ax1[1].set_xlabel('d (mm)')
    ax1[1].set_ylabel('P (d)')
    ax1[1].get_xaxis().get_major_formatter().set_useOffset(False)
    ax1[1].get_xaxis().get_major_formatter().set_scientific(False)
    ax1[1].set_xlim(data['d_lim'])
    plt.show(block=False)

    fig2, ax2 = plt.subplots(norder,nwaves,figsize=(12,10))
    axx = ax2.reshape(norder,nwaves)
    for i, w in enumerate(hists.keys()):
        for j, hist in enumerate(hists[w]):
            my_hist(axx[j,i], hist, bins=bins)
            axx[j,i].axvline(data['peaks'][w][j], color='k',zorder=15)
            axx[j,i].axvspan(data['peaks'][w][j]-data['pkerr'][w][j],data['peaks'][w][j]+data['pkerr'][w][j],color='gray',alpha=0.4,zorder=15)
            axx[j,i].set_ylabel('Order {0:d}'.format(data['orders'][w][j]))
            axx[j,i].get_xaxis().get_major_formatter().set_useOffset(False)
            axx[j,i].get_xaxis().get_major_formatter().set_scientific(False)
        axx[0,i].set_title('{0} nm'.format(w),fontsize=18)
        axx[-1,i].set_xlabel('R (px)')

    if saveit:
        fig0.savefig(join(fig_folder,'peaks.png'))
        fig1.savefig(join(fig_folder,'Ld_marginals.png'))
        fig2.savefig(join(fig_folder,'peak_histograms.png'))

    plt.show()


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
        parser.add_argument('--pkthresh', type=float, nargs='+', default=0.1,
                help='threshold for finding peak (can be a list corresponding to the various\
                peaks). Default is 0.1')
        parser.add_argument('--overwrite',action='store_true', help='allows you to overwrite\
                previously saved peak information')
        parser.add_argument('--pkerr', type=float, default=0.2, help='error for peak\
                locations in units of percent of binsize at peak. Should be between 1/6 and 1/2. Default is 1/5')
        parser.add_argument('--L_lim', '-L', type=float, nargs=2, default=[100.,110.], help='limit\
                for fitting L in units of mm. Default is 100-110.')
        parser.add_argument('--d_lim', '-d', type=float, nargs=2, default=[0.47,0.48], help='limit\
                for fitting d in units of mm. Default is 0.47-0.48.')
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
        
        if not isfile(fname) or args.overwrite:
            if isfile(org_fname):
                data = h5_2_dict(org_fname)
                print(data.keys())
                r = data['smooth_r']
                sig = data['smooth_sig']
                peaks, orders = get_pk_locations(r, sig, args.wavelengths, pkthresh=args.pkthresh)
                image_name = data['fname']
                if image_name[-3:].lower() == "nef":
                    image_data = get_data(image_name, color='b')
                    ny, nx = image_data.shape
                    x = np.arange(0, nx, 1)
                    y = np.arange(0, ny, 1)
                    x0, y0 = data['center']
                    xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
                    R = np.sqrt(xx**2 + yy**2)

                    R = R.flatten()
                    image_data = image_data.flatten()

                    idx_sort = np.argsort(R)

                    fig, ax = plt.subplots()
                    ax.plot(R[idx_sort], image_data[idx_sort])
                    for w in peaks:
                        for pk in peaks[w]:
                            ax.axvline(pk)
                    plt.show()




                print(peaks)
                print(orders)
                dict_2_h5(fname, {'r':r, 'sig':sig, 'peaks':peaks, 'orders':orders})
            else:
                raise ValueError('{0} does not exist!'.format(org_fname))
        else:
           a = h5_2_dict(fname)
           r = a['r']
           sig = a['sig']
           peaks = a['peaks']
           orders = a['orders']

        binsizes = np.diff(np.concatenate(([0],r)))
        pkerr = {}
        for i,w in enumerate(peaks.keys()):
            pkerrs = []
            for pk in peaks[w]:
                ix = np.abs(pk - r).argmin()
                pkerrs.append(args.pkerr * binsizes[ix])
            pkerr[w] = np.array(pkerrs)

        inputs = {'pkerr':pkerr, 'L_lim':args.L_lim, 
                'd_lim':args.d_lim, 'livepoints':args.livepoints}
        dict_2_h5(fname, inputs, append=True)

        if args.no_solve:
            solver_in = None
        else:
            solver_in = {'peaks':peaks, 'orders':orders, 'basename':basename, 'pkerr':pkerr,
                    'L_lim':args.L_lim, 'd_lim':args.d_lim, 'livepoints':args.livepoints, 'resume':resume}
    else:
        solver_in = None
    
    solver_in = Comm.bcast(solver_in, root=0)
    if solver_in is not None:
        ld_multinest_solver(solver_in['peaks'], solver_in['orders'], solver_in['basename'],
                solver_in['pkerr'], solver_in['L_lim'], solver_in['d_lim'], 
                livepoints=solver_in['livepoints'], resume=solver_in['resume'])

    if rank == 0:
        ld_check(args.folder,bins=args.bins)
