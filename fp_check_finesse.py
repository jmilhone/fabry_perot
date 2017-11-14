from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import plottingtools.core as ptools
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from os.path import join
import json
from histogram import my_hist
import fp_helpers
import argparse
import model2 as model
import time


rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

w = {'Th': 487.873302,
     'Ar': 487.98634}

mu = {'Th': 232.03806,
      'Ar': 39.948}


def make_hist_plots(F, A, rel, Ti, folder=None):
    fig, ax = plt.subplots()
    my_hist(ax, F, bins=20)
    ax.set_xlabel("F")
    ax.set_ylabel("P(F)")
    plt.show(block=False)

    fig1, ax1 = plt.subplots()
    my_hist(ax1, A, bins=20)
    ax1.ticklabel_format(style='sci', ax1is='x', scilimits=(0,0))
    ax1.set_xlabel("A (Counts)")
    ax1.set_ylabel("P(A)")
    plt.show(block=False)

    fig2, ax2 = plt.subplots()
    my_hist(ax2, rel, bins=20)
    ax2.set_xlabel("Rel. A")
    ax2.set_ylabel("P(Rel A.)")
    plt.show(block=False)

    fig3, ax3 = plt.subplots()
    my_hist(ax3, Ti, bins=20)
    ax3.set_xlabel("Ti Ar (eV)")
    ax3.set_ylabel("P(Ti Ar)")
    if folder:
        fig.savefig(join(folder, "F_marginal.pdf"))
        fig1.savefig(join(folder, "A_marginal.pdf"))
        fig2.savefig(join(folder, "Arel_marginal.pdf"))
        fig3.savefig(join(folder, "Ti_Ar_marginal.pdf"))
    plt.show(block=True)



def check_finesse_results(finesse_folder, Ld_folder, savefig=True):
    Lpost, dpost = fp_helpers.read_Ld_results(Ld_folder)
    Fpost, Amp_post, rel_post, Ti_post = fp_helpers.read_finesse_results(finesse_folder)

    folder = None
    if savefig:
        folder = join(finesse_folder, "Plots/")
        fp_helpers.make_directory(folder)

    make_hist_plots(Fpost, Amp_post, rel_post, Ti_post, folder=folder)
    LL = Lpost[::150]
    dd = dpost[::150]
    n = len(LL)

    FF = Fpost[::10]
    AA = Amp_post[::10]
    RR = rel_post[::10]
    Tii = Ti_post[::10]
    m = len(FF)

    data_fname = join(finesse_folder, "calib_data.json")
    with open(data_fname, 'r') as calibfile:
        calib_data = json.load(calibfile, parse_float=np.float64)

    r = np.array(calib_data['r'])
    s = np.array(calib_data['sig'])
    idx = calib_data['idx'][0]

    rr = r[idx]
    rr = np.linspace(rr.min()-50, rr.max()+50.0, 300)
    npts = len(rr)
    pred = np.zeros((npts, n, m))
    print npts, n, m

    Ti_Th = 0.025*1000.0/300.0
    for i in range(n):
        print i
        for j in range(m):
            pred[:, i, j] = model.forward4(rr, LL[i], dd[i], FF[j], [Ti_Th, Tii[j]],
                    [mu['Th'], mu['Ar']], [RR[j]*AA[j], AA[j]], [w['Th'], w['Ar']], V=0.0)

    pred = pred.reshape((npts, n*m))
    min_sig = np.min(pred, axis=1)
    max_sig = np.max(pred, axis=1)

    fig, ax = plt.subplots()
    ax.fill_between(rr, min_sig, max_sig, color='r', alpha=0.7)
    #ax.plot(r, s, 'b')
    ax.errorbar(r, s, yerr=0.03*s, fmt='.', color='b')
    ax.set_xlim(rr.min(), rr.max())
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    ax.set_xlabel("R (px)")
    ax.set_ylabel("Counts")
    if savefig:
        fig.savefig(join(folder, "Th_spectra_fit.pdf"))
    plt.show()
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve for the etalon spacing d and camera focal length L.")
    parser.add_argument("finesse_dir", metavar="Finesse directory", type=str, 
            help="Name of finesse directory located in Calibration/Ld_saves/")
    parser.add_argument("Ld_dir", metavar="Ld directory", type=str, 
            help="Name of the Ld directory with MultiNest results.")
    parser.add_argument("--order", "-O", action='store', type=int, 
            help="Order number to solve for", default=0)
    args = parser.parse_args()


    finesse_dir = join("Calibration/Finesse_saves/", args.finesse_dir, "")
    Ld_dir = join("Calibration/Ld_saves/", args.Ld_dir, "")
    calib_fname = join(finesse_dir, "calib_data.json")
    
    check_finesse_results(finesse_dir, Ld_dir)

