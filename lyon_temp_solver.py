import numpy as np
import argparse
from tools.file_io import h5_2_dict, dict_2_h5, prep_folder, read_Ld_results, read_match_finesse_results, read_lyon_temp_results
from tools.plotting import ringsum_click, my_hist, my_hist2d
from core.models import lyon_temp_forward
import matplotlib.pyplot as plt
from os.path import join, isfile, abspath
from mpi4py import MPI
import pymultinest
import subprocess
import itertools

def get_fit_region(r,s,error,plotit=True):
    
    edge_r,_ = ringsum_click(r,s,title='Click to the left and right of the fit fitting region')
    ix1 = np.abs(edge_r[0] - r).argmin()
    ix2 = np.abs(edge_r[1] - r).argmin()
    fit_ix = np.arange(ix1,ix2+0.5,1).astype(int)
    rr = r[fit_ix]
    ss = s[fit_ix]

    error = error * (ss.max()-ss.min()) * np.ones_like(ss) + error * ss

    if plotit:
        f,ax = plt.subplots()
        ax.plot(rr,ss)
        ax.fill_between(rr,  ss+error, ss-error, color='r', alpha=0.4)
        ax.set_xlim([rr.min()-0.1*(rr.max()-rr.min()),rr.max()+0.1*(rr.max()-rr.min())])
        #ax.plot(r,s-s.min())
        #pred = lyon_temp_forward(r, 205./0.004, 3.875, 4.5, 700, 0.15, -8000., 1.e8)
        #ax.plot(r,pred)
        plt.show()

    return fit_ix, error

def temp_solver(r, sig, sig_error, current, Lpost, dpost, Fpost, basename, livepoints=200, resume=True):
    T_lim = [0.01, 2.]
    V_lim = [-2., 2.]
    #A_lim = [sig.max(),1.e4*sig.max()]
    #logA_lim = [np.log10(x) for x in A_lim]
    #logO_lim = [-2., np.log10(10.*sig.min())]
    #A_lim = [0.5,2]

    ixs = range(len(Lpost))
    jxs = range(len(Fpost))

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(T_lim[1]-T_lim[0]) + T_lim[0]
        cube[1] = cube[1]*(V_lim[1]-V_lim[0]) + V_lim[0]
        #cube[2] = cube[2]*(A_lim[1]-A_lim[0]) + A_lim[0]
        #cube[2] = 10**(cube[2]*(logA_lim[1]-logA_lim[0]) + logA_lim[0])
        #cube[3] = 10**(cube[3]*(logO_lim[1]-logO_lim[0]) + logO_lim[0])

    def log_likelihood(cube, ndim, nparams):
        i = np.random.choice(ixs)
        j = np.random.choice(jxs)
        L = Lpost[i]
        d = dpost[i]
        F = Fpost[j]
        pred = lyon_temp_forward(r, L, d, F, current, cube[0], cube[1])#, E=Epost[j])
        pred /= pred.max()
        #pred *= cube[2]
        chisq = np.sum((sig - pred)**2 / sig_error**2)
        return -chisq / 2.0

    nparams = 2
    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=resume, verbose=True, sampling_efficiency='model', n_live_points=livepoints,
            outputfiles_basename=basename, max_modes=500)

def lyon_temp_check(folder, saveit=True, LdF_sub=20, T_sub=20, bins=30):
    
    if saveit:
        fig_folder = join(folder, 'temp_solver_plots')
        prep_folder(fig_folder)

    data = h5_2_dict(join(folder, 'input_temp_solver.h5'))
    
    Lpost, dpost = read_Ld_results(data['Ld_folder'])
    nLd = Lpost.size
    ix = np.random.choice(range(nLd), LdF_sub)
    LL = Lpost[ix]
    dd = dpost[ix]

    Fpost, _, _  = read_match_finesse_results(data['F_folder'])
    nF = Fpost.size
    ix = np.random.choice(range(nF), LdF_sub)
    FF = Fpost[ix]
    #EE = Epost[ix]

    Tpost, Vpost = read_lyon_temp_results(folder)
    nT = Tpost.size
    ix = np.random.choice(range(nT), T_sub)
    TT = Tpost[ix]
    VV = Vpost[ix]
    #AA = Apost[ix]

    rr = data['r'][data['fit_ix']]
    ss = data['sig'][data['fit_ix']]
    #ss -= ss.min()
    npts = len(rr)
    pred = np.zeros((npts, LdF_sub, T_sub))

    for i in range(LdF_sub):
        print '{0} of {1}'.format(i+1, LdF_sub)
        for j in range(T_sub):
            blah = lyon_temp_forward(rr, LL[i], dd[i], FF[i], data['current'],
                    TT[j], VV[j])#, E=EE[i])
            blah /= blah.max()
            #blah *= AA[j]
            pred[:, i, j] = blah

    pred = pred.reshape((npts, LdF_sub*T_sub))
    min_sig = np.min(pred, axis=1)
    max_sig = np.max(pred, axis=1)

    posts = [Tpost, Vpost]#, Apost]
    names = ['T (eV)', 'V (km/s)']#, 'A (counts)']

    fig, ax = plt.subplots(1, 2, figsize=(10,7))
    ax = ax.flatten()
    for i,p in enumerate(posts):
        my_hist(ax[i], p, bins=bins)
        ax[i].set_xlabel(names[i])
        ax[i].set_ylabel('P ({0})'.format(names[i].split(' ')[0]))
    fig.tight_layout()
    plt.show(block=False)

    fig1, ax1 = plt.subplots(figsize=(10,7))
    ax1.errorbar(rr, ss, yerr=data['error_per']*((ss.max()-ss.min())*np.ones_like(ss)+ss),fmt='.', color='b', zorder=1)
    ax1.fill_between(rr, min_sig, max_sig, color='r', alpha=0.55, zorder=10)
    ax1.set_xlabel('R (px)')
    ax1.get_xaxis().get_major_formatter().set_useOffset(False)
    ax1.get_xaxis().get_major_formatter().set_scientific(False)
    plt.show(block=False)

    #fig2, ax2 = plt.subplots(1, 2, figsize=(16,9))
    #ax2 = ax2.flatten()
    #names = list(itertools.combinations(names, 2))
    #for i, data in enumerate(itertools.combinations(posts,2)):
    #    my_hist2d(ax2[i], data[0], data[1], bins=bins)
    #    ax2[i].set_xlabel(names[i][1])
    #    ax2[i].set_ylabel(names[i][0])
    #fig2.tight_layout()

    if saveit:
        fig.savefig(join(fig_folder, 'marginals.png'))
        fig1.savefig(join(fig_folder, 'spectra.png'))
    #    fig2.savefig(join(fig_folder,'joint_marginals.png'))
    plt.show()

if __name__ == "__main__":
    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    if rank == 0:
        parser=argparse.ArgumentParser(description='Performs temp fit for lyon data, including Zeeman splitting')
        parser.add_argument('folder', type=str, help='folder with ringsum.h5 of image to use')
        parser.add_argument('ld_folder', type=str, help='ld calibration folder')
        parser.add_argument('match_finesse_folder', type=str, help='finesse calibration folder')
        parser.add_argument('current', type=float, help='current in coils for this shot, needed for Zeeman')
        parser.add_argument('--error', '-er', type=float, default=0.10, help='error percentage to use in fitting. default is 0.1')
        parser.add_argument('--overwrite', action='store_true', help='allows you to overwrite previously saved fitting region')
        parser.add_argument('--livepoints', type=int, default=100, help='number of livepoints to use in multinest. Default is 100.')
        parser.add_argument('--no_solve', action='store_true', help='only writes fit region information and, if it is there, plots output, skipping multinest call')
        parser.add_argument('--no_plot',action='store_true',help='supress plotting after multinest run, use for batch processing.')
        args = parser.parse_args()

        fname = abspath(join(args.folder, 'input_temp_solver.h5'))
        basename = abspath(join(args.folder, 'temp_'))
        org_fname = abspath(join(args.folder, 'ringsum.h5'))

        if args.overwrite:
            resume = False
        else:
            resume = True

        if not isfile(fname) or args.overwrite:
            if isfile(org_fname):
                data = h5_2_dict(org_fname)
                r = data['r']
                sig = data['sig']
                fit_ix, error = get_fit_region(r, sig, args.error)
                dict_2_h5(fname, {'r':r, 'sig':sig, 'fit_ix':fit_ix, 'error':error})
            else:
                raise ValueError('{0} does not exist!'.format(org_fname))
        else:
            a = h5_2_dict(fname)
            r = a['r']
            sig = a['sig']
            fit_ix = a['fit_ix']
            error = a['error']
        
        inputs = {'error_per':args.error, 'current':args.current, 'livepoints':args.livepoints,
                'Ld_folder':abspath(args.ld_folder), 'F_folder':abspath(args.match_finesse_folder)}
               
        dict_2_h5(fname, inputs, append=True)

        if args.no_solve:
            solver_in = None
        else:
            Lpost,dpost = read_Ld_results(abspath(args.ld_folder))
            Fpost,_,_ = read_match_finesse_results(abspath(args.match_finesse_folder))
            rr = r[fit_ix]
            ss = sig[fit_ix]
            #ss -= ss.min()
            solver_in = {'r':rr, 'sig':ss, 'error':error, 'current':args.current, #'Epost':Epost,
                    'Lpost':Lpost, 'dpost':dpost, 'Fpost':Fpost, 'resume':resume,
                    'basename':basename, 'livepoints':args.livepoints}
    else:
        solver_in = None

    solver_in = Comm.bcast(solver_in, root=0)
    if solver_in is not None:
        temp_solver(solver_in['r'], solver_in['sig'], solver_in['error'], solver_in['current'],
                solver_in['Lpost'], solver_in['dpost'], solver_in['Fpost'], #solver_in['Epost'],
                solver_in['basename'], livepoints=solver_in['livepoints'], resume=solver_in['resume'])
    
    if rank == 0:
        if not args.no_plot:
            lyon_temp_check(args.folder)

