from __future__ import division
from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
#import forward_model as fm
import model
import time
import plottingtools.core as ptools
import cPickle as pickle
import pymultinest
import json
import argparse
from os.path import join, isdir
from os import mkdir
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
### Globals     ###
### DO NOT EDIT ###
wTh = 487.873302
wAr = 487.98634

muAr = 40.0
muTh = 232.0
def calculate_peaks(L, d, w, norders=4):
    m = 2.e6 * d / w
    m0 = np.floor(m)
    #return [L * np.sqrt( m**2 / (m0 - j)**2 - 1.0) for j in range(norders)]
    return [np.sqrt(2.0) * L * np.sqrt(1.0 - (m0-j) / m) for j in range(norders)]


def read_L_d_results(savedir, basename="fp_Ld_"):

    analyzer = pymultinest.Analyzer(n_params=2, outputfiles_basename=join(savedir,basename))

    stats = analyzer.get_mode_stats()
    local_log_ev = [x['local log-evidence'] for x in stats['modes']]
    ix = np.argmax(local_log_ev)

    mode = stats['modes'][ix]
    mode_vals = mode['maximum a posterior']

    L = mode_vals[0]
    d = mode_vals[1]
    #print "L, d = ",L*.004, d
    return L, d



def retrieve_param_info(savedir, fname):
    f = join(savedir, fname)
    with open(f, 'r') as infile:
        param_info = json.load(infile)
    params = param_info['params']
    labels = param_info['labels']
    prob_labels = param_info['prob labels']
    return params, labels, prob_labels


def retrieve_input_data(savedir, fname, return_dict=False):
    f = join(savedir, fname)
    with open(f, 'rb') as infile:
        data = pickle.load(infile)
    ringsum = data['ringsum']
    r = data['binarr']
    idx = data['ind']
    if return_dict:
        return r, ringsum, idx, data
    else:
        return r, ringsum, idx

def check_finesse_solver(savedir, basename='fp_', param_fname="param_file.json", data_fname="fp_ringsum_params.p"):
    params, labels, prob_labels = retrieve_param_info(savedir, param_fname)
    r, ringsum, idx = retrieve_input_data(savedir, data_fname)

    analyzer = pymultinest.Analyzer(n_params=len(params), outputfiles_basename=join(savedir,basename))
    st = analyzer.get_stats()
    #print st['marginals']
    #print [x['1sigma'] for x in st['marginals']]
    stats = analyzer.get_mode_stats()
    log_ev = stats['global evidence']

    local_log_ev = [x['local log-evidence'] for x in stats['modes']]
    ix = np.argmax(local_log_ev)
    #mode = stats['modes'][0]
    
    mode = stats['modes'][ix]

    mode_vals = mode['maximum a posterior']

    F = mode_vals[0]
    Ti_Ar = mode_vals[1]
    Amp_Th = mode_vals[2]
    Amp_Ar = mode_vals[3]
    r0 = mode_vals[4]
    Ti_Th = 1000.0 * 0.025 / 300.0
    L, d = read_L_d_results("saves/Ld_test16")

    lin_data = model.forward4(r, L, d, F, [Ti_Th, Ti_Ar], [muTh, muAr], [Amp_Th, Amp_Ar],[wTh, wAr])
    #lin_data = model.forward4(r, L, d, F, [Ti_Th, Ti_Ar, Ti_Th], [muTh, muAr, muTh], [Amp_Th, Amp_Ar, .1*Amp_Ar],[wTh, wAr, 487.800942])
    lin_data *= np.exp(-(r/r0)**2)

    ringsum_sd = 0.01 * ringsum[idx] + 100.0 
    chisq = np.sum((lin_data[idx] - ringsum[idx])**2 / ringsum_sd**2)
    print "\nMaximum a Posterior Mode:"
    print "L = {} mm".format(L * .004)
    print "d = {} mm".format(d)
    print "Finesse = {}".format(F)
    print "Ti Ar = {} eV".format(Ti_Ar)
    print "Amp_Th = {} Counts".format(Amp_Th)
    print "Amp_Ar = {} Counts".format(Amp_Ar)
    print "R0 = {} pixels".format(r0)
    print "Chi Squared: {0:e}".format(chisq)
    print "Reduced Chi Squared: {0:e}".format(chisq / (len(idx)-len(params)))
    print "Log Evidence: {0:e}\n".format(log_ev)

    fig, ax = plt.subplots()
    #sigmas = ['1sigma', '2sigma', '3sigma']
    #alphas = [0.75, 0.5, 0.25]
    sigmas = ['3sigma']
    alphas = [1.0]
    rr = np.linspace(r.min(), r.max(), 10000)
    for alpha, sig in zip(alphas, sigmas):
        up = [x[sig][1] for x in st['marginals']]
        down = [x[sig][0] for x in st['marginals']]

        lin_up =model.forward4(rr, L, d, up[0], [Ti_Th, up[1]], [muTh, muAr], [up[2], up[3]],[wTh, wAr])  
        lin_up *= np.exp(-(rr/up[4])**2)

        lin_down =model.forward4(rr, L, d, down[0], [Ti_Th, down[1]], [muTh, muAr], [down[2], down[3]],[wTh, wAr])  
        lin_down *= np.exp(-(rr/down[4])**2)
        ax.fill_between(rr**2, lin_down, lin_up, facecolor='r', alpha=alpha)
        #ax.plot(r**2, lin_up)
        #ax.plot(r**2, lin_down)
    ax.plot(r**2, ringsum, 'b')
    ax.plot(r[idx]**2, ringsum[idx], 'g')
    ax.plot(r**2, lin_data, '--r')
    #ax.plot(r*0.004*350.0/150.0/25.4, ringsum, 'b')
    #ax.plot(r*0.004*350.0/150.0/25.4, lin_data, '--r')
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ptools.add_thick_box(ax, minor=False)
    ptools.add_labels(ax, r"R${}^2$ (px${}^2$)", "Counts")
    plt.tight_layout()
    plt.show()



def check_solver3(savedir, basename="fp_full_", param_fname="param_file.json", data_fname="fp_ringsum_params.p"):
    params, labels, prob_labels = retrieve_param_info(savedir, param_fname)
    r, ringsum, idx = retrieve_input_data(savedir, data_fname)

    analyzer = pymultinest.Analyzer(n_params=len(params), outputfiles_basename=join(savedir,basename))
    stats = analyzer.get_mode_stats()
    log_ev = stats['global evidence']

    local_log_ev = [x['local log-evidence'] for x in stats['modes']]
    ix = np.argmax(local_log_ev)
    #mode = stats['modes'][0]
    
    mode = stats['modes'][ix]

    mode_vals = mode['maximum a posterior']

    L = mode_vals[0]
    d = mode_vals[1]
    #print "messed with d"
    #d += 4.e-5
    F = mode_vals[2]
    Ti_Ar = mode_vals[3]
    Amp_Th = mode_vals[4]
    Amp_Ar = mode_vals[5]
    r0 = mode_vals[6]
    Ti_Th = 1000.0 * 0.025 / 300.0

    #lin_data = model.forward4(r, L, d, F, [Ti_Th, Ti_Ar], [muTh, muAr], [Amp_Th, Amp_Ar],[wTh, wAr])
    lin_data = model.forward4(r, L, d, F, [Ti_Th, Ti_Ar, Ti_Th], [muTh, muAr, muTh], [Amp_Th, Amp_Ar, .1*Amp_Ar],[wTh, wAr, 487.800942])
    lin_data *= np.exp(-(r/r0)**2)

    ringsum_sd = 0.01 * ringsum[idx] + 100.0 
    chisq = np.sum((lin_data[idx] - ringsum[idx])**2 / ringsum_sd**2)
    print "\nMaximum a Posterior Mode:"
    print "L = {} mm".format(L * .004)
    print "d = {} mm".format(d)
    print "Finesse = {}".format(F)
    print "Ti Ar = {} eV".format(Ti_Ar)
    print "Amp_Th = {} Counts".format(Amp_Th)
    print "Amp_Ar = {} Counts".format(Amp_Ar)
    print "R0 = {} pixels".format(r0)
    print "Chi Squared: {0:e}".format(chisq)
    print "Reduced Chi Squared: {0:e}".format(chisq / (len(idx)-len(params)))
    print "Log Evidence: {0:e}\n".format(log_ev)
    
    #w1 = 469.5038074  
    new_w = [469.5038074, 469.545425, 468.6194581, 468.335172, 468.925186, 468.925189, 469.062197, 
            467.3660927, 468.22319]
    for w1 in new_w:
        print  w1
        print calculate_peaks(L, d,w1 , norders=4)
        print [x**2 for x in calculate_peaks(L, d, w1, norders=4)]
        print ""
    #w2 =  469.545425 
    #print w2
    #print calculate_peaks(L, d,w2, norders=4)
    #print [x**2 for x in calculate_peaks(L, d, w2, norders=4)]
    fig, ax = plt.subplots()
    #ax.plot(r[idx]**2, ringsum[idx], 'b')
    ax.plot(r**2, ringsum, 'b')
    ax.plot(r**2, lin_data, '--r')
    #ax.plot(r*0.004*350.0/150.0/25.4, ringsum, 'b')
    #ax.plot(r*0.004*350.0/150.0/25.4, lin_data, '--r')
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ptools.add_thick_box(ax, minor=False)
    ptools.add_labels(ax, r"R${}^2$ (px${}^2$)", "Counts")
    plt.tight_layout()
    plt.show()


def check_solver2(savedir, basename="fp_newfull_", param_fname="param_file.json", data_fname="fp_ringsum_params.p"):
    params, labels, prob_labels = retrieve_param_info(savedir, param_fname)
    r, ringsum, idx, data = retrieve_input_data(savedir, data_fname, return_dict=True)

    analyzer = pymultinest.Analyzer(n_params=len(params), outputfiles_basename=join(savedir,basename))
    stats = analyzer.get_mode_stats()
    local_log_ev = [x['local log-evidence'] for x in stats['modes']]
    ix = np.argmax(local_log_ev)
    log_ev = stats['global evidence']

    mode = stats['modes'][0]
    mode = stats['modes'][ix]

    mode_vals = mode['maximum a posterior']

    L = mode_vals[0]
    d = mode_vals[1]
    F = mode_vals[2]
    Ti_Ar = mode_vals[3]
    rel_amp = data['rel_amp']
    Amp_Ar = mode_vals[4]
    r0 = data['rscale']
    Ti_Th = 1000.0 * 0.025 / 300.0

    lin_data = model.forward2(r, L, d, F, [Ti_Th, Ti_Ar], [muTh, muAr], [Amp_Ar / rel_amp, Amp_Ar],[wTh, wAr])
    lin_data *= np.exp(-(r/r0)**2)

    ringsum_sd = 0.01 * ringsum[idx] + 100.0 
    chisq = np.sum((lin_data[idx] - ringsum[idx])**2 / ringsum_sd**2)
    print "\nMaximum a Posterior Mode:"
    print "L = {} mm".format(L * .004)
    print "d = {} mm".format(d)
    print "Finesse = {}".format(F)
    print "Ti Ar = {} eV".format(Ti_Ar)
    print "Amp_Ar = {} Counts".format(Amp_Ar)
    print "Chi Squared: {0:e}".format(chisq)
    print "Reduced Chi Squared: {0:e}".format(chisq / (len(idx)-len(params)))
    print "Log Evidence: {0:e}\n".format(log_ev)

    fig, ax = plt.subplots()
    ax.plot(r**2, ringsum, 'b')
    ax.plot(r**2, lin_data, '--r')
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ptools.add_thick_box(ax, minor=False)
    ptools.add_labels(ax, r"R${}^2$ (px${}^2$)", "Counts")
    plt.tight_layout()
    plt.show()


def check_full_solver(savedir, basename="fp_full_", param_fname="param_file.json", data_fname="fp_ringsum_params.p"):
    params, labels, prob_labels = retrieve_param_info(savedir, param_fname)
    r, ringsum, idx = retrieve_input_data(savedir, data_fname)

    analyzer = pymultinest.Analyzer(n_params=len(params), outputfiles_basename=join(savedir,basename))
    stats = analyzer.get_mode_stats()
    local_log_ev = [x['local log-evidence'] for x in stats['modes']]
    ix = np.argmax(local_log_ev)
    #mode = stats['modes'][0]
    mode = stats['modes'][ix]

    log_ev = stats['global evidence']
    #mode_vals = mode['maximum a posterior']
    mode_vals = mode['maximum a posterior']

    L = mode_vals[0]
    d = mode_vals[1]
    F = mode_vals[2]
    Ti_Th = mode_vals[3]
    Ti_Ar = mode_vals[4]
    Amp_Th = mode_vals[5]
    Amp_Ar = mode_vals[6]
    r0 = mode_vals[7]

    lin_data = model.forward2(r, L, d, F, [Ti_Th, Ti_Ar], [muTh, muAr], [Amp_Th, Amp_Ar],[wTh, wAr])
    lin_data *= np.exp(-(r/r0)**2)

    ringsum_sd = 0.01 * ringsum[idx] + 100.0 
    chisq = np.sum((lin_data[idx] - ringsum[idx])**2 / ringsum_sd**2)

    print "\nMaximum a Posterior Mode:"
    print "L = {} mm".format(L * .004)
    print "d = {} mm".format(d)
    print "Finesse = {}".format(F)
    print "Ti Th = {} K".format(Ti_Th * 1000.0 / .025)
    print "Ti Ar = {} eV".format(Ti_Ar)
    print "Amp_Th = {} Counts".format(Amp_Th)
    print "Amp_Ar = {} Counts".format(Amp_Ar)
    print "R0 = {} pixels".format(r0)
    print "Chi Squared: {0:e}".format(chisq)
    print "Reduced Chi Squared: {0:e}".format(chisq / (len(idx)-len(params)))
    print "Log Evidence: {0:e}\n".format(log_ev)

    fig, ax = plt.subplots()
    ax.plot(r**2, ringsum, 'b')
    ax.plot(r**2, lin_data, '--r')
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ptools.add_thick_box(ax, minor=False)
    ptools.add_labels(ax, r"R${}^2$ (px${}^2$)", "Counts")
    plt.tight_layout()
    plt.show()

def plot_marginals(savedir, basename="fp_full_",param_fname="param_file.json", save=False):
    params, labels, prob_labels = retrieve_param_info(savedir, param_fname)
    f = basename + "post_equal_weights.dat"
    f = join(savedir, f)
    with open(f, 'r') as infile:
        post = np.loadtxt(infile, ndmin=2)
    plt.plot(post[:,1])
    plt.show()
    folder = join(savedir, "Plots/")
    if save:
        if not isdir(folder):
            try:
                mkdir(folder)
            except OSError, e:
                print "Error making directory"
                raise e

    for idx, param in enumerate(params):
        fig, ax = plt.subplots()
        if params[idx].lower() == 'l':
            hist, bins = np.histogram(post[:, idx]*.004, density=True, bins=40)
        else:
            hist, bins = np.histogram(post[:, idx], density=True, bins=40)
        bw = bins[1]-bins[0]
        ax.bar(bins[0:-1], hist*bw, width=bw)
        ptools.add_labels(ax, labels[idx], prob_labels[idx])
        ptools.add_thick_box(ax, minor=False)

        plt.tight_layout()
        if save:
            fname = join(folder, "{}_marginal.pdf".format(param))
            fig.savefig(fname)
            plt.close(fig)
        else:
            plt.show()

    for i, param1 in enumerate(params):
        for j, param2 in enumerate(params[0:i]):
            fig, ax = plt.subplots()
            hist, xx, yy = np.histogram2d(post[:, j], post[:, i], bins=100, normed=True)
            dx = xx[1]-xx[0]
            dy = yy[1]-yy[0]
            im = ax.contourf(yy[0:-1], xx[0:-1], hist*dx*dy)
            ax_divider = make_axes_locatable(ax)
            cax = ax_divider.append_axes("right", size="7%", pad="2%")
            cb = plt.colorbar(im, cax=cax)
            cax.tick_params(labelsize=18)
            ptools.add_thick_box(ax, minor=False)
            ptools.add_labels(ax, labels[i], labels[j])
            plt.tight_layout()
            if save:
                fname = join(folder, "{}_{}_marginal.pdf".format(param2, param1))
                fig.savefig(fname)
                plt.close(fig)
            else:
                plt.show()

if __name__ == "__main__":
    #check_solver3("saves/solver3_run0", basename="fp_full_", param_fname="param_file.json", data_fname="fp_ringsum_params.p")
    #check_full_solver("saves/full_solver_run17", basename="fp_full_", param_fname="param_file.json", data_fname="fp_ringsum_params.p")

    parser = argparse.ArgumentParser()
    parser.add_argument("--savedir", "-s", action='store', dest='savedir', type=str, 
            #default="saves/finesse_solver8", help="Directory where MultiNest save files are located")
            default="saves/finesse_solver10/", help="Directory where MultiNest save files are located")
    parser.add_argument("--basename", "-b", action='store', dest='base_name', type=str, 
            #default="fp_", help="Basename for MultiNest save file")
            default="fp_", help="Basename for MultiNest save file")
    parser.add_argument("--paramname", "-p", action='store', dest='param_name', type=str,
            default="param_file.json", help="Name for file containing parameters used in MultiNest solution")
    parser.add_argument("--dataname", "-d", action='store', dest='data_name', type=str,
            default="fp_ringsum_params.p", help="Name for file containing ring sum data for MultiNest solver")
    parser.add_argument("--nparams", "-n", action='store', dest='nparams', type=int, default=7, help="Number of parameters used in MultiNest") 
    args = parser.parse_args()

    check_finesse_solver(args.savedir, basename=args.base_name, param_fname=args.param_name, data_fname=args.data_name)
    #plot_marginals("saves/full_solver_run17", basename="fp_full_", save=False)
    #if args.nparams == 7:
    #    check_solver3(args.savedir, basename=args.base_name, param_fname=args.param_name, data_fname=args.data_name)
    ##elif args.nparams == 8:
    ##    check_full_solver(args.savedir, basename=args.base_name, param_fname=args.param_name, data_fname=args.data_name)
    #elif args.nparams == 5:
    #    check_solver2(args.savedir, basename=args.base_name, param_fname=args.param_name, data_fname=args.data_name)



    #plot_marginals(args.savedir, basename=args.base_name,param_fname=args.param_name, save=False)
