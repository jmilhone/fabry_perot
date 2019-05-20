from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np
from os.path import join, abspath
from fabry.tools import images, plotting, file_io
from fabry.core import fitting, ringsum
import skimage.io as io
import argparse
from tqdm import tqdm


def interpolate_point(pt1, pt2, thres):
    x1, y1 = pt1
    x2, y2 = pt2
    slope = (y2-y1)/(x2-x1)
    offset = y2 - slope*x2

    point = (thres - offset) / slope 
    return point

def half_max_finder(x, y, max_val_loc, halfmax):
    loc = np.abs(x-max_val_loc).argmin()
    iR = np.min(np.where(y[loc:] < halfmax)) + loc
    iL = np.max(np.where(y[:loc] < halfmax))

    L = interpolate_point((x[iL], y[iL]), (x[iL+1], y[iL+1]), halfmax)
    R = interpolate_point((x[iR-1], y[iR-1]), (x[iR], y[iR]), halfmax)

    return L, R


def quick_gaussian_fit(x, y):
    p0 = np.polyfit(x, np.log(y), 2)
    sigma = np.sqrt(-0.5 / p0[0])
    mu = -p0[1] / 2 / p0[0]
    lnA = p0[2] - p0[1]**2 / 4 / p0[0]
    Amp = np.exp(lnA)

    return Amp, mu, sigma

def determine_fit_indicies(x, y, pk_range = [150, 250]):

    left = np.abs(x-pk_range[0]).argmin()
    right = np.abs(x-pk_range[1]).argmin()+1
    sl = slice(left, right)
    ipk = np.argmax(y[sl]) + left
    pk = x[ipk]**2
    indicies = fitting.determine_fit_range(x**2, y, pk, thres=0.2)

    return indicies



def main(fname, binsize=0.2):

    data = io.imread(fname, plugin='tifffile')
    nimages, nrow, ncol = data.shape
    idx = 8
    xguess, yguess = plotting.center_plot(data[idx, :, :])

    centers = list() 
    centers.append(ringsum.locate_center(data[idx, :, :], xguess=xguess, 
        yguess=yguess, block_center=False, binsize=0.1, plotit=False, printit=True))
    #nimages=5
    #for i in range(1, nimages):
    #    centers.append(ringsum.locate_center(data[i, :, :], 
    #        xguess=centers[i-1][0], yguess=centers[i-1][1], block_center=False,
    #        binsize=0.1, plotit=False, printit=False))
    #    print(centers[-1])


    #for center in centers:
    #    print(center)

    r_list = list()
    sig_list =list()
    sd_list = list()
    for i in tqdm(range(nimages)):
        x0, y0 = centers[0]
        r, sig, sd = ringsum.ringsum(data[i,:,:],x0,y0, use_weighted=False, 
                quadrants=False, binsize=binsize, remove_hot_pixels=True) 
        r_list.append(r)
        sig_list.append(sig)
        sd_list.append(sd)
    print('subtracting the first image')
    sig_list = [x - sig_list[0] for x in sig_list]
    sd_list = [np.sqrt(x**2 + sd_list[0]**2) for x in sd_list]

    peaks = []
    fwhm = []
    for i in range(nimages):
        r = r_list[i]
        sig = sig_list[i]
        if sig.max() > 10:
            j = determine_fit_indicies(r, sig, pk_range=[150, 250])
            amp, mu, sigma = quick_gaussian_fit(r[j]**2, sig[j])
            mu = np.sqrt(mu)
            peaks.append(mu)
            LL, RR = half_max_finder(r**2, sig, mu**2, amp/2.0)
            LL = np.sqrt(LL)
            RR = np.sqrt(RR)
            fwhm.append(RR-LL)
        else:
            fwhm.append(np.nan)
            peaks.append(np.nan)

    fig, ax = plt.subplots(2, sharex=True)
    ax[0].plot(peaks, 'o')
    ax[1].plot(fwhm, 'o')
    plt.show()

    #n = len(r_list[0])
    #r = r_list[0]
    #sig = np.zeros((len(r_list), n))
    #sig_sd = np.zeros_like(sig)
    #for idx, (ssig, ssd) in enumerate(zip(sig_list, sd_list)):
    #    sig[idx, :] = ssig / ssig.max()
    #    sig_sd[idx, :] = ssd
    markers = ['-', '-.', '--']
    fig, ax = plt.subplots()
    for idx, (r, sig, sd) in enumerate(zip(r_list, sig_list, sd_list)):
        n, m = divmod(idx, 10)
        ax.errorbar(r, sig, yerr=sd, linestyle=markers[n], label=str(idx))
        print(idx, sig.max())
        # if sig.max() > 10:
        #     ax.plot(r, sig/sig.max(), markers[n], label=str(idx))
    ax.legend()
    #ax.errorbar(r, np.mean(sig, axis=0), yerr=np.std(sig, axis=0))
    #ax.plot(r, np.std(sig, axis=0) / np.mean(sig, axis=0))
    #print([len(x) for x in r_list])
    ax.set_xlim(160, 260)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='stores image ringsum into Data folder')
    parser.add_argument('fname',type=str,help='NEF filename of image you wish to process')
    args = parser.parse_args()
    main(args.fname)

