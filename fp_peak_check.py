from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import plottingtools.core as ptools
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from os.path import join
import json
from histogram import my_hist
from fp_helpers import peak_calculator, make_directory
import argparse
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

wavelengths = {'Th': 487.873302,
               'Ar': 487.98634}


def check_peaks(folder, peak_fname, savefig=True):

    if savefig:
        fig_folder = join(folder, "Plots/")
        make_directory(fig_folder)

    with open(peak_fname, 'r') as peak_file:
        peak_data = json.load(peak_file, parse_float=np.float64)

    orders = {}
    nAr = len(peak_data['Ar'])
    nTh = len(peak_data['Th'])
    n = max((nAr, nTh))

    orders['Ar'] = np.linspace(0.0, nAr-1.0, nAr)
    orders['Th'] = np.linspace(0.0, nTh-1.0, nTh)

    fname = join(folder, "fp_Ld_post_equal_weights.dat")
    post = np.loadtxt(fname, ndmin=2)
    npts = post.shape[0]


    Ar_peaks = np.zeros((nAr, npts))
    Th_peaks = np.zeros((nTh, npts))


    for idx, order in enumerate(orders['Ar']):
        Ar_peaks[idx, :] = peak_calculator(post[:, 0], post[: ,1], wavelengths['Ar'], order)

    for idx, order in enumerate(orders['Th']):
        Th_peaks[idx, :] = peak_calculator(post[:, 0], post[: ,1], wavelengths['Th'], order)

    fig, ax = plt.subplots()
    my_hist(ax, post[:, 0]*0.004, bins=100)
    ax.set_xlabel("L (mm)")
    ax.set_ylabel("P(L)")
    plt.show(block=False)


    fig1, ax1 = plt.subplots()
    my_hist(ax1, post[:, 1], bins=100)
    ax1.set_xlabel("d (mm)")
    ax1.set_ylabel("P(d)")
    plt.show(block=False)


    fig2, ax2 = plt.subplots(n, 2)
    for i in range(nTh):
        my_hist(ax2[i][0], Th_peaks[i, :], bins=100)
        ax2[i][0].axvline(peak_data['Th'][i], color='k')
        ax2[i][0].set_ylabel("Order {0:d}".format(i))

    for i in range(nAr):
        my_hist(ax2[i][1], Ar_peaks[i, :], bins=100)
        ax2[i][1].axvline(peak_data['Ar'][i], color='k')
    ax2[0][0].set_title("Th I 487.87 nm")
    ax2[0][1].set_title("Ar II 487.98 nm")
    ax2[-1][0].set_xlabel("R (px)")
    ax2[-1][1].set_xlabel("R (px)")
    fig.tight_layout()

    if savefig:
        fig.savefig(join(fig_folder, "L_marginal.pdf"))
        fig1.savefig(join(fig_folder, "d_marginal.pdf"))
        fig2.savefig(join(fig_folder, "Peak_histograms.pdf"))
    plt.show()


    #w = [487.800942, 487.649495, 488.08585, 488.09142, 488.120436, 488.185311]
    #for ww in w:
    #    fig, ax = plt.subplots(4)
    #    for i in range(4):
    #        pk = peak_calculator(post[:, 0], post[:, 1], ww, i)
    #        my_hist(ax[i], pk, bins=100)
    #    ax[0].set_title(str(ww))
    #    ax[-1].set_xlabel("R (px)")
    #    plt.show(block=False)
    #plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solve for the etalon spacing d and camera focal length L.")
    parser.add_argument("Ld_dir", metavar="Ld directory", type=str, 
            help="Name of Ld directory located in Calibration/Ld_saves/")
    args = parser.parse_args()

    pkname = "calibration_peaks.json"

    Ld_dir = join("Calibration/Ld_saves/", args.Ld_dir, "")
    pkname = join(Ld_dir, pkname)

    check_peaks(Ld_dir, pkname)

