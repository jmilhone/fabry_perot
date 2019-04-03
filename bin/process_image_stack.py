from __future__ import print_function, division
import matplotlib.pyplot as plt
import numpy as np
from os.path import join, abspath
from fabry.tools import images, plotting, file_io
from fabry.core import fitting, ringsum
import skimage.io as io
import argparse
from tqdm import tqdm


def main(fname, binsize=0.1):

    data = io.imread(fname, plugin='tifffile')
    nimages, nrow, ncol = data.shape
    xguess, yguess = plotting.center_plot(data[0, :, :])

    centers = list() 
    centers.append(ringsum.locate_center(data[0, :, :], xguess=xguess, 
        yguess=yguess, block_center=False, binsize=0.1, plotit=False, printit=True))
    #nimages=10
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

    #n = len(r_list[0])
    #r = r_list[0]
    #sig = np.zeros((len(r_list), n))
    #sig_sd = np.zeros_like(sig)
    #for idx, (ssig, ssd) in enumerate(zip(sig_list, sd_list)):
    #    sig[idx, :] = ssig / ssig.max()
    #    sig_sd[idx, :] = ssd
    fig, ax = plt.subplots()
    for r, sig, sd in zip(r_list, sig_list, sd_list):
        ax.errorbar(r, sig, yerr=sd)
    #    ax.plot(r, sig/sig.max())
    #ax.errorbar(r, np.mean(sig, axis=0), yerr=np.std(sig, axis=0))
    #ax.plot(r, np.std(sig, axis=0) / np.mean(sig, axis=0))
    #print([len(x) for x in r_list])
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='stores image ringsum into Data folder')
    parser.add_argument('fname',type=str,help='NEF filename of image you wish to process')
    args = parser.parse_args()
    main(args.fname)

