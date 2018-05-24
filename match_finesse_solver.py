import numpy as np
import argparse
from tools.file_io import h5_2_dict, dict_2_h5, prep_folder, read_Ld_results, read_match_finesse_results
from tools.plotting import ringsum_click, peak_plot, my_hist, my_hist2d
import matplotlib.pyplot as plt
from core.models import match_finesse_forward,lyon_temp_forward
from os.path import join,isfile,abspath
from mpi4py import MPI
import pymultinest
import subprocess
import itertools

def get_finesse_region(r,s,error,plotit=True):

    edge_r,_ = ringsum_click(r,s,title='Click to the left and right of the finesse fitting region')
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
        plt.show()

    return fit_ix, error

def temp_match_solver(r, sig, sig_error, temp, temp_sigma, Lpost, dpost, basename, F_lim=(1.,18.), V_lim = [-2., 2.], livepoints=500, w0=487.98634, mu=39.948, errtemp=True,resume=True, arb_error=0.0):

    #A_lim = [0.5, 2]
    T_lim = [0.025, temp+10.*temp_sigma]
    ixs = range(len(Lpost))

    if errtemp:
        E_lim = [0.05, 0.5]
        nparams = 4
        def log_prior(cube, ndim, nparams):
            cube[0] = cube[0]*(F_lim[1]-F_lim[0]) + F_lim[0]
            cube[1] = cube[1]*(V_lim[1]-V_lim[0]) + V_lim[0]
            cube[2] = cube[2]*(T_lim[1]-T_lim[0]) + T_lim[0]
            cube[3] = cube[3]*(E_lim[1]-E_lim[0]) + E_lim[0]
            #cube[4] = cube[4]*(A_lim[1]-A_lim[0]) + A_lim[0]

        def log_likelihood(cube, ndim, nparams):
            j = np.random.choice(ixs)
            L = Lpost[j]
            d = dpost[j]
            pred = match_finesse_forward(r, L, d, cube[0], cube[2], cube[1], errtemp=cube[3], w0=w0, mu=mu)
            pred /= pred.max()
            #pred *= cube[4]
            chisq = np.sum((sig - pred)**2 / sig_error**2) + ((cube[2] - temp)**2 / temp_sigma**2)
            return -chisq / 2.0
    else:
        nparams = 3
        def log_prior(cube, ndim, nparams):
            cube[0] = cube[0]*(F_lim[1]-F_lim[0]) + F_lim[0]
            cube[1] = cube[1]*(V_lim[1]-V_lim[0]) + V_lim[0]
            cube[2] = cube[2]*(T_lim[1]-T_lim[0]) + T_lim[0]
            #cube[3] = cube[3]*(A_lim[1]-A_lim[0]) + A_lim[0]

        def log_likelihood(cube, ndim, nparams):
            j = np.random.choice(ixs)
            L = Lpost[j]
            d = dpost[j]
            #pred = match_finesse_forward(r, L, d, cube[0], cube[2], cube[1], errtemp=None, w0=w0, mu=mu)
            pred = lyon_temp_forward(r, L, d, cube[0], 50, cube[2], cube[1]) 
            pred /= pred.max()
            #pred *= cube[3]
            chisq = np.sum((sig - pred)**2 / sig_error**2) + ((cube[2] - temp)**2 / temp_sigma**2)
            return -chisq / 2.0

    pymultinest.run(log_likelihood, log_prior, nparams, importance_nested_sampling=False,
            resume=resume, verbose=True, sampling_efficiency='model', n_live_points=livepoints,
            outputfiles_basename=basename, max_modes=500)

def match_finesse_check(folder, saveit=True, Ld_sub=20, F_sub=20, bins=30):

    if saveit:
        fig_folder = join(folder, 'match_finesse_plots')
        prep_folder(fig_folder)
    
    data = h5_2_dict(join(folder,'input_match_finesse.h5'))
    Lpost, dpost = read_Ld_results(data['Ld_folder'])
    nLd = Lpost.size
    ix = np.random.choice(range(nLd), Ld_sub)
    LL = Lpost[ix]
    dd = dpost[ix]
    
    
    if data['errtemp']:
        Fpost, Vpost, Tpost, Epost = read_match_finesse_results(folder, errtemp=True)
        nF = Fpost.size
        ix = np.random.choice(range(nF), F_sub)
        FF = Fpost[ix]
        VV = Vpost[ix]
        #AA = Apost[ix]
        #OO = Opost[ix]
        TT = Tpost[ix]
        EE = Epost[ix]
    else:
        #Fpost, Vpost, Apost, Opost, Tpost = read_match_finesse_results(folder)
        Fpost, Vpost, Tpost = read_match_finesse_results(folder)
        nF = Fpost.size
        ix = np.random.choice(range(nF), F_sub)
        FF = Fpost[ix]
        VV = Vpost[ix]
        #AA = Apost[ix]
        #OO = Opost[ix]
        TT = Tpost[ix]

    if 'fit_ix' not in data.keys():
        rr = data['r']
        ss = data['sig']
    else:
        rr = data['r'][data['fit_ix']]
        ss = data['sig'][data['fit_ix']]
    #ss -= ss.min()
    npts = len(rr)
    pred = np.zeros((npts, Ld_sub, F_sub))

    if data['errtemp']:
        for i in range(Ld_sub):
            print '{0} of {1}'.format(i+1, Ld_sub)
            for j in range(F_sub):
                blah = match_finesse_forward(rr, LL[i], dd[i], 
                          FF[j], TT[j], VV[j], errtemp=EE[j],w0=data['wavelength'],mu=data['mu'])# + OO[j]
                blah /= blah.max()
                #blah *= AA[j]
                pred[:, i, j] = blah
    else:
        for i in range(Ld_sub):
            print '{0} of {1}'.format(i+1, Ld_sub)
            for j in range(F_sub):
                #blah = match_finesse_forward(rr, LL[i], dd[i], 
                #          FF[j], TT[j], VV[j],errtemp=None,w0=data['wavelength'],mu=data['mu'])# + OO[j]
                blah = lyon_temp_forward(rr, LL[i], dd[i], FF[j], 50, TT[j], VV[j])
                blah /= blah.max()
                #blah *= AA[j]
                pred[:,i,j] = blah

    pred = pred.reshape((npts, Ld_sub*F_sub))        
    min_sig = np.min(pred, axis=1)
    max_sig = np.max(pred, axis=1)

    if data['errtemp']:
        posts = [Fpost, Vpost, Tpost, Epost]
        names = ['F', 'V (km/s)', r'T$_{LIF}$ (eV)', r'T$_{extra}$ (eV)']
    else:
        posts = [Fpost, Vpost, Tpost]
        names = ['F', 'V (km/s)', r'T$_{LIF}$ (eV)']

    fig, ax = plt.subplots(2,2,figsize=(10,7))
    ax = ax.flatten()
    for i,p in enumerate(posts):
        my_hist(ax[i], p, bins=bins)
        ax[i].set_xlabel(names[i])
        ax[i].set_ylabel('P ({0})'.format(names[i].split(' ')[0]))
        if i == 2:
            _,binz = np.histogram(p, density=True, bins=bins)
            bw = binz[1]-binz[0]
            ti_x = np.linspace(data['temp']-10*data['temp_sigma'],data['temp']+10*data['temp_sigma'],1000)
            ti_prior = np.exp(-0.5*(ti_x - data['temp'])**2/(data['temp_sigma'])**2)/(data['temp_sigma']*np.sqrt(2.*np.pi))
            xlim = ax[i].get_xlim()
            ax[i].plot(ti_x,ti_prior*bw,'r',label='prior')
            ax[i].set_xlim(xlim) 
            ax[i].legend(frameon=False)

    fig.tight_layout()
    plt.show(block=False)

    fig1,ax1 = plt.subplots(figsize=(12,10))
    ax1.errorbar(rr,ss,yerr=data['error_per']*((ss.max()-ss.min())*np.ones_like(ss)+ss),fmt='.',color='b',zorder=1)
    ax1.fill_between(rr, min_sig, max_sig, color='r', alpha=0.7, zorder=10)
    ax1.set_xlabel('R (px)')
    ax1.get_xaxis().get_major_formatter().set_useOffset(False)
    ax1.get_xaxis().get_major_formatter().set_scientific(False)
    plt.show(block=False)

    if data['errtemp']:
        fig2,ax2 = plt.subplots(2,3,figsize=(16,9))
    else:
        fig2,ax2 = plt.subplots(2,2,figsize=(16,6))
    ax2 = ax2.flatten()
    names = list(itertools.combinations(names,2))
    for i,data in enumerate(itertools.combinations(posts,2)):
        my_hist2d(ax2[i],data[0],data[1],bins=bins)
        ax2[i].set_xlabel(names[i][1])
        ax2[i].set_ylabel(names[i][0])
    fig2.tight_layout()

    if saveit:
        fig.savefig(join(fig_folder,'marginals.png'))
        fig1.savefig(join(fig_folder,'spectra_fit.png'))
        fig2.savefig(join(fig_folder,'joint_marginals.png'))
    plt.show()
    
if __name__ == "__main__":
    Comm = MPI.COMM_WORLD
    rank = Comm.Get_rank()
    if rank == 0:
        parser = argparse.ArgumentParser(description='Performs finesse calibration matching with a known temperature')
        parser.add_argument('folder', type=str, help='folder of finesse image to use, this is also where multinest data is saved')
        parser.add_argument('ld_folder', type=str, help='folder containing Ld calibration multinest output')
        parser.add_argument('--temp','-t',type=float,default=0.15,help='known temperature to use for finesse fitting. Default is 0.15 eV (which is the case for Lyon measurements with 700W of power and 50A of coil current according to LIF')
        parser.add_argument('--temp_sigma','-ts',type=float,default=0.05,help='known temperature standard deviation. Default is 0.05(0.015) eV')
        parser.add_argument('--wavelength','-w0',type=float,default=487.98634,help='wavelength of line in nm, default is Ar II 487.98634')
        parser.add_argument('--mass','-mu',type=float, default=39.948,help='atomic mass of atom/ion producing the line (for doppler broadening calc) in amu, default is Argon=>39.948 others: Helium=>4.002, Thorium 232.03806')
        parser.add_argument('--error','-er',type=float,default=0.10,help='error percentage to use in fitting. default is 0.10')
        parser.add_argument('--arb_error','-aer',type=float,default=0.0,help='error percentage to arbitrarily add to the chisq calc. default is 0')
        parser.add_argument('--extra_gauss', '-eg', action='store_true', help='flag to add extra broadening mechinism to forward model')
        parser.add_argument('--overwrite',action='store_true',help='allows you to overwrite previously saved finesse region')
        parser.add_argument('--F_lim', '-F', type=float, nargs=2, default=[1.,10.], help='bounds for finesse fit. Default is 1-10')
        parser.add_argument('--V_lim', '-V', type=float, nargs=2, default=[-2.,2.], help='bounds for velocity fit (km/s). Default is -2->2')
        parser.add_argument('--livepoints', type=int, default=500, help='number of livepoints to use in multinest. Default is 500.')
        parser.add_argument('--no_solve', action='store_true', help='only writes finesse region information, skipping multinest solve.')
        args=parser.parse_args()

        fname = abspath(join(args.folder, 'input_match_finesse.h5'))
        basename = abspath(join(args.folder, 'F_'))
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
                fit_ix, error = get_finesse_region(r, sig, args.error)
                dict_2_h5(fname, {'r':r, 'sig':sig, 'fit_ix':fit_ix, 'error':error})
            else:
                raise ValueError('{0} does not exist!'.format(org_fname))
        else:
            a = h5_2_dict(fname)
            r = a['r']
            sig = a['sig']
            fit_ix = a['fit_ix']
            error = a['error']
        
        inputs = {'error_per':args.error, 'F_lim':args.F_lim, 'temp':args.temp, 'wavelength':args.wavelength,'mu':args.mass, 'arb_error':args.arb_error,
                'livepoints':args.livepoints, 'Ld_folder':abspath(args.ld_folder), 'errtemp':args.extra_gauss, 'V_lim':args.V_lim, 'temp_sigma':args.temp_sigma}
               
        dict_2_h5(fname, inputs, append=True)

        if args.no_solve:
            solver_in = None
        else:
            Lpost,dpost = read_Ld_results(abspath(args.ld_folder))
            rr = r[fit_ix]
            ss = sig[fit_ix]
            #ss -= ss.min()
            solver_in = {'r':rr,'sig':ss,'basename':basename, 'error':error, 'resume':resume, 'arb_error':args.arb_error,
                    'F_lim':args.F_lim, 'livepoints':args.livepoints, 'V_lim': args.V_lim, 'temp_sigma': args.temp_sigma,
                    'temp':args.temp, 'w0':args.wavelength, 'mu':args.mass, 'Lpost':Lpost, 'dpost':dpost, 'errtemp':args.extra_gauss}
    else:
        solver_in = None

    solver_in = Comm.bcast(solver_in, root=0)
    if solver_in is not None:
        temp_match_solver(solver_in['r'], solver_in['sig'], solver_in['error'], 
                solver_in['temp'], solver_in['temp_sigma'], solver_in['Lpost'], solver_in['dpost'], 
                solver_in['basename'], F_lim=solver_in['F_lim'], V_lim=solver_in['V_lim'], errtemp=solver_in['errtemp'],
                livepoints=solver_in['livepoints'], w0=solver_in['w0'], mu=solver_in['mu'], resume=solver_in['resume'], arb_error=solver_in['arb_error'])
    
    if rank == 0:
        match_finesse_check(args.folder)
