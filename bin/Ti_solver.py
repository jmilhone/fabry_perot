from __future__ import print_function, division
import sys
sys.path.append("../")
import numpy as np
import argparse
from fabry.tools.file_io import h5_2_dict, dict_2_h5, prep_folder, read_Ld_results, read_finesse_results
from fabry.tools.plotting import ringsum_click, my_hist, tableau20_colors
import matplotlib.pyplot as plt
from fabry.core.models import forward_model
from os.path import join,isfile,abspath
from mpi4py import MPI
import pymultinest
import multiprocessing as mp
import random 

def get_fitting_region(r,s,error,plotit=True):

    edge_r,_ = ringsum_click(r,s,title='Click to the left and right of the fitting region')
    ix1 = np.abs(edge_r[0] - r).argmin()
    ix2 = np.abs(edge_r[1] - r).argmin()
    fit_ix = np.arange(ix1,ix2+0.5,1).astype(int)
    rr = r[fit_ix]
    ss = s[fit_ix]

    # error = error * np.mean(ss) * np.ones_like(ss) + error * ss

    if plotit:
        f,ax = plt.subplots()
        ax.plot(rr,ss)
        ax.fill_between(rr,  ss+error[fit_ix], ss-error[fit_ix], color='r', alpha=0.4)
        ax.set_xlim([rr.min()-0.1*(rr.max()-rr.min()),rr.max()+0.1*(rr.max()-rr.min())])
        plt.show()

    return fit_ix, None


def Ti_solver(r, sig, sig_error, Ld_dir, finesse_dir, Ti_lim, V_lim, A_lim, basename, resume=False, 
        w0=487.98634, mu=39.948, livepoints=1000):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0]*(Ti_lim[1] - Ti_lim[0]) + Ti_lim[0]
        cube[1] = 1000.0*(cube[1]*(V_lim[1] - V_lim[0]) + V_lim[0]) # need to convert to m/s from km/s
        cube[2] = cube[2]*(A_lim[1] - A_lim[0]) + A_lim[0]

    def log_likelihood(cube, ndim, nparams):
        #i = np.random.choice(nL)
        #j = np.random.choice(nF)

        #L = Lpost[i]
        #d = dpost[i]
        #F = Fpost[j]

        L = 0.380173301412519577E+05
        d = 0.883628502371783142E+00
        F = 21.0
        vals = forward_model(r, L, d, F, w0, mu, cube[2], cube[0], cube[1], sm_ang=False, nlambda=1024)
        chisq = np.sum((vals - sig)**2 / sig_error**2)

        return -chisq/2.0
    r = r[::3]
    sig = sig[::3]
    sig_error = sig_error[::3]
    nparams = 3

    Lpost, dpost = read_Ld_results(Ld_dir)
    Fpost, _, _, _ = read_finesse_results(finesse_dir)

    nL = len(Lpost)
    nF = len(Fpost)

    #L = 0.380173301412519577E+05
    #d = 0.883628502371783142E+00
    #F = 21.0
    #npts = 200
    #test_sig = np.zeros((npts, len(r)))
    #for i in xrange(npts):
    #    cube = [random.random() for _ in xrange(3)]
    #    log_prior(cube, None, None)
    #    test_sig[i, :] = forward_model(r, L, d, F, w0, mu, cube[2], cube[0], cube[1], sm_ang=False, nlambda=1024)
    #fig, ax = plt.subplots()
    #for i in xrange(npts):
    #    ax.plot(r, test_sig[i, :], 'C0')
    #ax.errorbar(r, sig, yerr=sig_error, color='C1', ecolor='C2')
    #plt.show()



    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=resume, verbose=True, sampling_efficiency='model', n_live_points=livepoints,
            outputfiles_basename=basename)


def check_solution(folder, Ld_dir, finesse_dir, recover=False, w0=487.98634, mu=39.948):
    print("I'm here!")
    data = h5_2_dict(join(folder, 'input_plasma_data.h5'))

    ix = data['fit_ix']
    r = data['r']
    sig = data['sig']

    Lpost, dpost = read_Ld_results(Ld_dir)
    Fpost, _, _, _ = read_finesse_results(finesse_dir)
    nL = len(Lpost)
    nF = len(Fpost)
    i = np.random.choice(nL)
    j = np.random.choice(nF)

    L = Lpost[i]
    d = dpost[i]
    F = Fpost[j]

    Lstep = 100
    Fstep = 100
    Tistep = 100
    analyzer = pymultinest.Analyzer(3, outputfiles_basename=join(folder, "Ti_"))
    modes = analyzer.get_mode_stats()
    Ti, V, A = modes['modes'][0]['mean']
    print(modes['modes'][0]['sigma'])
    print(Ti, V, A)
    post = analyzer.get_equal_weighted_posterior()

    Tipost = post[::, 0]
    Vpost = post[::, 1]
    Apost = post[::, 2]

    if recover:
        try:
            post_dict = h5_2_dict(join(folder, "Ti_solver_model_post.h5"))
            sig_post = post_dict["signal post"]
        except:
            print("Can't recover Ti solver model posterior.  Calculating from scratch.")
            sig_post = calculate_signal_post(r[ix], Lpost[::Lstep], dpost[::Lstep], Fpost[::Fstep],
                Tipost[::Tistep], Vpost[::Tistep], Apost[::Tistep], w0, mu, nprocs=32)
            dict_2_h5(join(folder, "Ti_solver_model_post.h5"), {"signal post": sig_post})
    else:
        sig_post = calculate_signal_post(r[ix], Lpost[::Lstep], dpost[::Lstep], Fpost[::Fstep],
             Tipost[::Tistep], Vpost[::Tistep], Apost[::Tistep], w0, mu, nprocs=32)
        dict_2_h5(join(folder, "Ti_solver_model_post.h5"), {"signal post": sig_post})


    sig_mean = np.mean(sig_post, axis=1)
    percentiles = calculate_percentile_ranges(sig_post)
    #vals = forward_model(r[ix], L, d, F, w0, mu, A, Ti, V, sm_ang=False, nlambda=1024)

    fig, ax = plt.subplots(figsize=(3.5, 3.5/1.61))
    # ax.plot(r[ix], sig[ix], 'C1', alpha=0.5, label='Data')
    ax.plot(r[ix], sig[ix]/100.0, 'C1', alpha=0.5, label='Data')
    #ax.plot(r[ix], vals, 'r')
    alphas = [0.8, 0.5, 0.2]
    keys = [68, 95, 99]
    for alpha, per in zip(alphas, keys):
        if per == 99:
            # ax.fill_between(r[ix], percentiles[per][0], percentiles[per][1], color='C3', alpha=alpha, label='Fit')
            ax.fill_between(r[ix], percentiles[per][0]/100.0, percentiles[per][1]/100.0, color='C3', alpha=alpha, label='Fit')
        else:
            # ax.fill_between(r[ix], percentiles[per][0], percentiles[per][1], color='C3', alpha=alpha)
            ax.fill_between(r[ix], percentiles[per][0]/100.0, percentiles[per][1]/100.0, color='C3', alpha=alpha)
    fig.legend(frameon=False, fontsize=8, loc='upper right', bbox_to_anchor=(0.5, 0.5))
    ax.set_xlabel("R (px)", fontsize=8, labelpad=-1)
    ax.set_ylabel("Counts (Hundreds)", fontsize=8, labelpad=-1)
    ax.tick_params(labelsize=8)
    #fig.tight_layout()
    fig.savefig(join(folder, "Ti_Ar_fit.png"), dpi=400)
    plt.show(block=False)

    axis_labels = ["Ti (eV)", "V (m/s)","A (Counts)"]
    ylabels = ["P(Ti)", "P(V)", "P(A)"]
    # fig, ax = plt.subplots(3, figsize=(6, 15))
    fig, ax = plt.subplots(2, figsize=(3.5, 2*3.5/1.61))
    # for n in range(3):
    for n in range(2):
        my_hist(ax[n], post[:, n])
        ax[n].set_xlabel(axis_labels[n], fontsize=8, labelpad=-1)
        ax[n].set_ylabel(ylabels[n], fontsize=8, labelpad=-1)
        ax[n].tick_params(labelsize=8)
    #fig.tight_layout()
    fig.savefig(join(folder, "Ti_solver_histograms.png"), dpi=400)
    plt.show()


def calculate_percentile_ranges(data):
    levels = [68, 95, 99]
    lower_levels = [50.0 - x/2.0 for x in levels]
    upper_levels = [50.0 + x/2.0 for x in levels]

    percentiles = {}
    for level, LL, UL in zip(levels, lower_levels, upper_levels):
        lower = np.percentile(data, LL, axis=1)
        upper = np.percentile(data, UL, axis=1)
        percentiles[level] = (lower, upper)

    return percentiles


def calculate_signal_post(r, Lpost, dpost, Fpost, Tipost, Vpost, Apost, w0, mu, nprocs=16):
    Lp, Fp, Tip = np.meshgrid(Lpost, Fpost, Tipost)
    dp, _, Vp = np.meshgrid(dpost, Fpost, Vpost)
    _, _, Ap = np.meshgrid(Lpost, Fpost, Apost)

    Lp = Lp.flatten()
    dp = dp.flatten()
    Fp = Fp.flatten()
    Tip = Tip.flatten()
    Vp = Vp.flatten()
    Ap = Ap.flatten()
    print(len(Lp))

    Lps = np.array_split(Lp, nprocs)
    dps = np.array_split(dp, nprocs)
    Fps = np.array_split(Fp, nprocs)
    Tips = np.array_split(Tip, nprocs)
    Vps = np.array_split(Vp, nprocs)
    Aps = np.array_split(Ap, nprocs)

    procs = []
    sigs = {}
    out = mp.Queue()
    labels = ["{0}".format(proc) for proc in range(nprocs)]

    for k in range(nprocs):
        p = mp.Process(target=_forward_model, args=(r, Lps[k], dps[k], Fps[k], Aps[k], Tips[k], Vps[k], w0, mu),
                kwargs={'out': out, 'label': labels[k]})
        procs.append(p)
        p.start()

    for i in range(nprocs):
        tup = out.get()
        sigs[tup[0]] = tup[1]

    for p in procs:
        p.join()

    sig_post = None
    for k in labels:
        if sig_post is None:
            sig_post = sigs[k]
        else:
            sig_post = np.hstack((sig_post, sigs[k]))
    return sig_post


def _forward_model(r, Lpost, dpost, Fpost, Apost, Tipost, Vpost, w0, mu, out=None, label=None):
    n = len(Lpost)
    sig = np.zeros((len(r), n))
    for i in range(n):
        sig[:, i] = forward_model(r, Lpost[i], dpost[i], Fpost[i], w0, mu, Apost[i], 
                Tipost[i], Vpost[i], sm_ang=False, nlambda=1024) 

    if out and label:
        out.put((label, sig))
    else:
        return sig



if __name__ == "__main__":
    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    if rank == 0:
        parser = argparse.ArgumentParser(description='Performs a Ti and V solver in Ar.')

        parser.add_argument('folder', type=str,
                help='folder of finesse image to use, this is also where multinest data is saved')

        parser.add_argument('ld_folder', type=str, help='folder containing Ld calibration multinest output')

        parser.add_argument('finesse_folder', type=str,
                help='folder containing finesse calibration multinest output')

        parser.add_argument('--Ti_lim', '-Ti', type=float, nargs=2, default=[0.005,5.0],
                help='bounds for Ti fit. Default is 0.005-5 eV')

        parser.add_argument('--V_lim', type=float, nargs=2, default=[-2.0, 2.0],
                help='bounds for V fit.  Default is -2 to +2 km/s')

        parser.add_argument('--overwrite',action='store_true',
                help='allows you to overwrite previously saved finesse region')

        parser.add_argument('--no_solve', action='store_true',
                help='only writes finesse region information, skipping multinest solve.')

        parser.add_argument('--error','-er',type=float,default=0.05,
                help='error percentage to use in fitting. default is 0.05')
        
        parser.add_argument('--livepoints', type=int, default=1000,
                help='number of livepoints to use in multinest. Default is 500.')

        parser.add_argument('--recover', action='store_true',
                help=("Recover finesse solver model posterior written to an h5 file because "
                      "calculation takes a long time"))
        args = parser.parse_args()

        fname = abspath(join(args.folder, "input_plasma_data.h5"))
        basename = abspath(join(args.folder, 'Ti_'))
        org_fname = abspath(join(args.folder, 'ringsum.h5'))

        if not isfile(fname) or args.overwrite:
            if isfile(org_fname):
                data = h5_2_dict(org_fname)
                r = data['r']
                sig = data['sig']
                error = data['sig_sd']
                fit_ix, _ = get_fitting_region(r, sig, error)
                maxval = np.max(sig[fit_ix])
                amp_lim = [0.5*maxval, 5.0*maxval]
                dict_2_h5(fname, {'r':r, 'sig':sig, 'fit_ix':fit_ix, 'error':error})
            else:
                raise ValueError('{0} does not exist!'.format(org_fname))
        else:
            a = h5_2_dict(fname)
            r = a['r']
            sig = a['sig']
            fit_ix = a['fit_ix']
            error = a['error']
            maxval = np.max(sig[fit_ix])
            amp_lim = [0.5*maxval, 5.0*maxval]

        inputs = {'error_per': args.error, 'Ti_lim': args.Ti_lim, 'V_lim': args.V_lim, 
                'livepoints': args.livepoints, 'Ld_folder':abspath(args.ld_folder), 
                'finesse_folder':abspath(args.finesse_folder), "A_lim": amp_lim}

        dict_2_h5(fname, inputs, append=True)

        if args.no_solve:
            solver_in = None
        else:
            rr = r[fit_ix]
            ss = sig[fit_ix]
            ss_sd = error[fit_ix]
            solver_in = {'r':rr,'sig':ss,'basename':basename, 'error':ss_sd, 
                    "Ld_dir": abspath(args.ld_folder), 'finesse_dir': abspath(args.finesse_folder),
                    "Ti_lim": args.Ti_lim, "V_lim": args.V_lim, "A_lim": amp_lim, 'w0': 487.98634,
                    "mu": 39.948, 'livepoints': args.livepoints}
    else:
        solver_in = None

    solver_in = Comm.bcast(solver_in, root=0)
    if solver_in is not None:
        print('im cheating!')
        Ti_solver(solver_in['r'], solver_in['sig']-150.0, solver_in['error'], solver_in['Ld_dir'],
                solver_in['finesse_dir'], solver_in['Ti_lim'], solver_in['V_lim'], solver_in['A_lim'],
                solver_in['basename'], resume=False, w0=solver_in['w0'], mu=solver_in['mu'], 
                livepoints=solver_in['livepoints'],)


    if rank == 0:
        check_solution(args.folder, abspath(args.ld_folder), abspath(args.finesse_folder), recover=args.recover)













