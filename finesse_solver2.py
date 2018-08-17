from __future__ import print_function, division
import numpy as np
import argparse
from tools.file_io import h5_2_dict, dict_2_h5, prep_folder, read_Ld_results, read_finesse_results
from tools.plotting import ringsum_click, my_hist, tableau20_colors
import matplotlib.pyplot as plt
from core.models import forward_model
from os.path import join,isfile,abspath
from mpi4py import MPI
import pymultinest
import multiprocessing as mp
from core.models import peak_calculator

def get_finesse_region(r,s,error,plotit=True):

    edge_r,_ = ringsum_click(r,s,title='Click to the left and right of the finesse fitting region')
    ix1 = np.abs(edge_r[0] - r).argmin()
    ix2 = np.abs(edge_r[1] - r).argmin()
    fit_ix = np.arange(ix1,ix2+0.5,1).astype(int)
    rr = r[fit_ix]
    ss = s[fit_ix]

    error = error * np.mean(ss) * np.ones_like(ss) + error * ss

    if plotit:
        f,ax = plt.subplots()
        ax.plot(rr,ss)
        ax.fill_between(rr,  ss+error, ss-error, color='r', alpha=0.4)
        ax.set_xlim([rr.min()-0.1*(rr.max()-rr.min()),rr.max()+0.1*(rr.max()-rr.min())])
        plt.show()

    return fit_ix, error


def finesse_solver(r, sig, sig_error, Lpost, dpost, basename, F_lim, A_lim, 
        Arel_lim, Ar_Ti_lim, offset_lim, w0, mu, livepoints=500, resume=True):
        #Arel_lim, Ar_Ti_lim,  w0, mu, livepoints=500):

    def log_prior(cube, ndim, nparams):
        cube[0] = cube[0] * (F_lim[1] - F_lim[0]) + F_lim[0]
        cube[1] = cube[1] * (A_lim[1] - A_lim[0]) + A_lim[0]
        cube[2] = cube[2] * (Arel_lim[1] - Arel_lim[0]) + Arel_lim[0]
        cube[3] = cube[3] * (Ar_Ti_lim[1] - Ar_Ti_lim[0]) + Ar_Ti_lim[0]
        cube[4] = 10**(cube[4]*(offset_log_lim[1] - offset_log_lim[0]) + offset_log_lim[0])

    def log_likelihood(cube, ndim, nparams):
        # forward_model(r, L, d, F, w0, mu, amp, temp, v, nlambda=1024, sm_ang=True)
        i = np.random.choice(npost)
        L = Lpost[i]
        d = dpost[i]
        vals = forward_model(r, L, d, cube[0], w0, mu,
                [cube[2]*cube[1], cube[1]], [0.025*1000.0/300.0, cube[3]],
                [0.0, 0.0], sm_ang=False, nlambda=1024)
        vals += cube[4]
        chisq = np.sum((vals - sig)**2 / sig_error**2)
        return -chisq/2.0

    offset_log_lim = [np.log10(x) for x in offset_lim]
    nparams = 5
    npost = len(Lpost)
    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=resume, verbose=True, sampling_efficiency='q', n_live_points=livepoints,
            outputfiles_basename=basename, max_modes=500)

def check_finesse(folder, recover=False):
    data = h5_2_dict(join(folder, 'input_finesse.h5'))
    ix = data['fit_ix'] 
    Lpost, dpost = read_Ld_results(data['Ld_folder'])

    i = np.random.choice(len(Lpost))
    L = Lpost[i]
    d = dpost[i]
    r = data['r']
    sig = data['sig']

    analyzer = pymultinest.Analyzer(5, outputfiles_basename=join(folder, "finesse_"))
    modes = analyzer.get_mode_stats()
    F, A, Arel, Ti, offset = modes['modes'][0]['mean']
    #F, A, Arel, Ti = modes['modes'][0]['mean']
    post = analyzer.get_equal_weighted_posterior()
    w0 = [487.873302, 487.98634]
    mu = [232.03806, 39.948]

    Lpost = Lpost[::30]
    dpost = dpost[::30]

    Fpost = post[::30, 0]
    Apost = post[::30, 1]
    Arelpost = post[::30, 2]
    Tipost = post[::30, 3]
    offsetpost = post[::30, 4]
    
    offset_mean = np.mean(post[:, 4])
    offset_sd = np.std(post[:, 4])
    print('F: {0:f} +/- {1:f}'.format(np.mean(post[:, 0]), np.std(post[:, 0])))
    print('A: {0:e} +/- {1:e}'.format(np.mean(post[:, 1]), np.std(post[:, 1])))
    print('Arel: {0:f} +/- {1:f}'.format(np.mean(post[:, 2]), np.std(post[:, 2])))
    print('Ti Ar: {0:f} +/- {1:f}'.format(np.mean(post[:, 3]), np.std(post[:, 3])))
    print('offset: {0:e} +/- {1:e}'.format(offset_mean, offset_sd))

    if recover:
        try:
            post_dict = h5_2_dict(join(folder, "finesse_solver_model_post.h5"))
            sig_post = post_dict["signal post"] 
        except IOError, e:
            print("Can't recover finesse solver q posterior. Calculating from scratch.")
            sig_post = calculate_signal_post(r[ix], Lpost, dpost, Fpost, Apost, Arelpost, Tipost, offsetpost, w0, mu)
            dict_2_h5(join(folder, "finesse_solver_model_post.h5"), {'signal post': sig_post})
    else:
        sig_post = calculate_signal_post(r[ix], Lpost, dpost, Fpost, Apost, Arelpost, Tipost, offsetpost, w0, mu)
        dict_2_h5(join(folder, "finesse_solver_model_post.h5"), {'signal post': sig_post})

    sig_mean = np.mean(sig_post, axis=1)
    sig_std = np.std(sig_post, axis=1)
    percentiles = calculate_percentile_ranges(sig_post)

    fig, ax = plt.subplots()
    ax.plot(r[ix], sig[ix], 'C0', alpha=0.7)
    #ax.plot(r[ix], sig_mean, 'C3')
    alphas = [0.8, 0.5, 0.2]
    #alphas = [0.9, 0.7, 0.5]
    keys = [68, 95, 99]
    for alpha, per in zip(alphas, keys):
        ax.fill_between(r[ix], percentiles[per][0], percentiles[per][1], color='C3', alpha=alpha)
    ax.set_xlabel("R (px)", fontsize=18)
    ax.set_ylabel("Counts", fontsize=18)
    ax.tick_params(labelsize=16)
    fig.tight_layout()
    plt.show(block=False)

    axis_labels = ["F", "A (Counts)", "Arel", "Ti (eV)", "Offset (Counts)"]
    ylabels = ["P(F)", "P(A)", "P(Arel)", "P(Ti)"]
    # fig1, ax1 = plt.subplots(3, 2)
    fig1, ax1 = plt.subplots(2, 2)
    for n in range(4):
        i, j = divmod(n, 2)
        my_hist(ax1[i][j], post[:, n])
        ax1[i][j].set_xlabel(axis_labels[n], fontsize=16)
        ax1[i][j].set_ylabel(ylabels[n], fontsize=16)
        ax1[i][j].tick_params(labelsize=14)
    fig1.tight_layout()
    plt.show(block=False)

    fig, ax = plt.subplots()
    ax.plot(r**2, sig-np.min(sig), 'C0', alpha=0.7)
    #ax.axhspan(offset_mean-offset_sd, offset_mean+offset_sd, color='C2', alpha=0.7)
    #ax.plot(r[ix], sig_mean, 'C3')
    alphas = [0.8, 0.5, 0.2]
    #alphas = [0.9, 0.7, 0.5]
    keys = [68, 95, 99]
    for alpha, per in zip(alphas, keys):
        ax.fill_between(r[ix]**2, percentiles[per][0], percentiles[per][1], color='C3', alpha=alpha)
    ax.set_xlabel("R^2 (px)", fontsize=18)
    ax.set_ylabel("Counts", fontsize=18)
    ax.tick_params(labelsize=16)

    orders = 1
    # w = [487.62247, 487.649495, 488.08585, 488.09142, 488.185311, 487.800942, 488.15961]
    w = [488.185311, 487.800942, 488.15961]
    amps = [0.087, 0.13, 0.065]

    L, d = read_Ld_results("Data/2018_04_11/Argon1/")
    idx = len(L)
    idx = np.random.choice(idx)
    Lmean = L[idx]
    dmean = d[idx]

    Fmean = np.mean(post[:, 0])
    Amean = np.mean(post[:, 1])
    Arelmean = np.mean(post[:, 2])
    Timean = np.mean(post[:, 3])

    print(Lmean, dmean, Fmean, Amean)
    amps = [x * Amean for x in amps]
    mus = [232.0 for x in range(len(amps))] 
    temps = [1000.0 * 0.025 / 300.0 for x in range(len(amps))]
    vs = [0.0 for x in range(len(amps))]

    w += [487.873302, 487.98634]
    amps += [Amean * Arelmean, Amean]
    mus += [232.0, 40.0]
    temps += [temps[0], Timean]
    vs += [0.0, 0.0]
    rr = np.linspace(0.0, 1000.0, 1001)
    newoffset = forward_model(rr, Lmean, dmean, Fmean, w, mus, amps, temps, vs, sm_ang=False)
    for i, w0 in enumerate(w):
        for order in range(orders):
            pks  = peak_calculator(L, d, w0, order)
            ax.axvline(np.mean(pks)**2, label=str(w0), color="C{0:d}".format(i+1))
    ax.legend()
    #ax.plot(rr**2, newoffset, 'k')
    fig.tight_layout()


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

def calculate_signal_post(r, Lpost, dpost, Fpost, Apost, Arelpost, Tipost, offsetpost, w0, mu, nprocs=32):
    # This is ugly....
    procs = []
    sigs = {}
    out = mp.Queue()
    labels = ["{0}".format(proc) for proc in range(nprocs)]
    nL = len(Lpost)
    nF = len(Fpost)
    Lp = Lpost[:, np.newaxis] * np.ones((nL, nF))
    dp = dpost[:, np.newaxis] * np.ones((nL, nF))

    Fp = np.ones((nL, nF)) * Fpost[np.newaxis, :]
    Ap = np.ones((nL, nF)) * Apost[np.newaxis, :]
    Arelp = np.ones((nL, nF)) * Arelpost[np.newaxis, :]
    Tip = np.ones((nL, nF)) * Tipost[np.newaxis, :]
    offsetp = np.ones((nL, nF)) * offsetpost[np.newaxis, :]

    Lp = Lp.flatten()
    dp = dp.flatten()
    Fp = Fp.flatten()
    Ap = Ap.flatten()
    Arelp = Arelp.flatten()
    Tip = Tip.flatten()
    offsetp = offsetp.flatten()

    Lps = np.array_split(Lp, nprocs)
    dps = np.array_split(dp, nprocs)
    Fps = np.array_split(Fp, nprocs)
    Aps = np.array_split(Ap, nprocs)
    Arelps = np.array_split(Arelp, nprocs)
    Tips = np.array_split(Tip, nprocs)
    offsetps = np.array_split(offsetp, nprocs)

    for k in range(nprocs):
        p = mp.Process(target=_forward_model, args=(r, Lps[k],  dps[k], Fps[k], Aps[k], Arelps[k], Tips[k], offsetps[k], w0, mu),
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


def _forward_model(r, Lpost, dpost, Fpost, Apost, Arelpost, Tipost, offsetpost, w0, mu, out=None, label=None):
    n = len(Lpost)
    sig = np.zeros((len(r), n))
    for i in range(n):
        sig[:, i] = forward_model(r, Lpost[i], dpost[i], Fpost[i], w0, mu, [Apost[i]*Arelpost[i], Apost[i]], 
                [0.025*1000.0/300.0, Tipost[i]], [0.0, 0.0], sm_ang=False, nlambda=1024) 
        sig[:, i] += offsetpost[i]
    if out and label:
        out.put((label, sig))
    else:
        return sig


if __name__ == "__main__":
    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    if rank == 0:
        parser = argparse.ArgumentParser(description='Performs finesse calibration in Ar.')
        parser.add_argument('folder', type=str, 
                help='folder of finesse image to use, this is also where multinest data is saved')
        parser.add_argument('ld_folder', type=str, help='folder containing Ld calibration multinest output')
        parser.add_argument('--F_lim', '-F', type=float, nargs=2, default=[15.,25.], 
                help='bounds for finesse fit. Default is 15-25')
        parser.add_argument('--Arel_lim', type=float, nargs=2, default=[0.1, 0.7], 
                help='bounds for relative amplitude of Th I to Ar II.  Default is 0.1-0.7')
        parser.add_argument('--Ti_Ar_lim', '-T', type=float, nargs=2, default=[0.025, 1.0], 
                help='bounds for Ar Ti fit.  Default is 0.025-1.0 eV')
        parser.add_argument('--overwrite',action='store_true',
                help='allows you to overwrite previously saved finesse region')
        parser.add_argument('--no_solve', action='store_true', 
                help='only writes finesse region information, skipping multinest solve.')
        parser.add_argument('--error','-er',type=float,default=0.05,
                help='error percentage to use in fitting. default is 0.05')
        parser.add_argument('--livepoints', type=int, default=500, 
                help='number of livepoints to use in multinest. Default is 500.')
        parser.add_argument('--recover', action='store_true', 
                help=("Recover finesse solver q posterior written to an h5 file because "
                      "calculation takes a long time"))
        args = parser.parse_args()

        fname = abspath(join(args.folder, "input_finesse.h5"))
        basename = abspath(join(args.folder, 'finesse_'))
        org_fname = abspath(join(args.folder, 'ringsum.h5'))

        if not isfile(fname) or args.overwrite:
            if isfile(org_fname):
                data = h5_2_dict(org_fname)
                r = data['r']
                sig = data['sig']
                fit_ix, error = get_finesse_region(r, sig, args.error)
                maxval = np.max(sig[fit_ix])
                amp_lim = [0.1*maxval, 10.0*maxval]
                offset_lim = [0.1, np.abs(np.min(sig))]
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
            offset_lim = [0.1, np.abs(np.min(sig))]
            amp_lim = [0.1*maxval, 10.0*maxval]
        print(offset_lim)
        inputs = {'error_per':args.error, 'F_lim':args.F_lim, 'Ti_lim':args.Ti_Ar_lim, 'Arel_lim': args.Arel_lim,
                'livepoints':args.livepoints, 'Ld_folder':abspath(args.ld_folder), 'Amp_lim': amp_lim,
                'offset_lim': offset_lim,
                 }

        dict_2_h5(fname, inputs, append=True)

        if args.no_solve:
            solver_in = None
        else:
            Lpost, dpost = read_Ld_results(abspath(args.ld_folder))
            rr = r[fit_ix]
            ss = sig[fit_ix]
            solver_in = {'r':rr,'sig':ss,'basename':basename, 'error':error,
                    'F_lim':args.F_lim, 'livepoints':args.livepoints, 'Ti_lim': args.Ti_Ar_lim, 
                    'Arel_lim': args.Arel_lim, 'Amp_lim':amp_lim, 'w0':[487.873302, 487.98634], 
                    'mu':[232.03806, 39.948], 'Lpost':Lpost, 'dpost':dpost, 'offset_lim': offset_lim}
    else:
        solver_in = None

    solver_in = Comm.bcast(solver_in, root=0)
    if solver_in is not None:
        finesse_solver(solver_in['r'], solver_in['sig'], solver_in['error'], solver_in['Lpost'], 
                solver_in['dpost'], solver_in['basename'], solver_in['F_lim'], solver_in['Amp_lim'], 
                solver_in['Arel_lim'], solver_in['Ti_lim'], solver_in['offset_lim'], solver_in['w0'], 
                solver_in['mu'], livepoints=solver_in['livepoints'], resume=True)


    if rank == 0:
        check_finesse(args.folder, recover=args.recover)





