from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
import time

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'




if __name__ == "__main__":
    binsize = 0.1
    w = 487.98634
    n_dlambda = 1536
    delta_lambda = 0.0001

    folder = "Images"
    fname = join(folder, "0015676_000.nef")
    bg_fname = join(folder, "0015676_001.nef")

    data = im.get_image_data(fname, bg_fname, color='b')

    L = 37562.1081426  # pixels
    d = .879999589322  # mm

    x0, y0 = (3067.7, 2037.6)
    x0, y0 = rs.locate_center(data, x0, y0, binsize=binsize, plotit=False)

    binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)

    peaks = fitting.peak_and_fit(binarr, sigarr, thres=0.55, plotit=False)
    r1 = peaks[0]
    r2 = peaks[1]
    m1 = 2.0 * d * 1.0e6 * L / np.sqrt(L ** 2 + r1 ** 2) / w
    m2 = 2.0 * d * 1.0e6 * L / np.sqrt(L ** 2 + r2 ** 2) / w

    print m1, m2

    ny, nx = data.shape
    x = np.arange(1, nx + 1, 1)
    y = np.arange(1, ny + 1, 1)
    xx, yy = np.meshgrid(1. * x - x0, 1. * y - y0)
    R = np.sqrt(xx ** 2 + yy ** 2)

    # Flatten data for simplicity
    R = R.flatten()
    data = data.flatten()

    lambda_arr = np.arange(-n_dlambda / 2, n_dlambda / 2, 1) * delta_lambda + w
    lambda_min = lambda_arr.min()
    lambda_max = lambda_arr.max()
    npeaks = len(peaks)

    ringsums = rs.proper_ringsum(R, data, m1, L, d, peaks, lambda_min,
                                 lambda_max, delta_lambda, ndl=n_dlambda)
    wpeaks = []
    for idx in range(npeaks):
        wpeaks += fitting.peak_and_fit(lambda_arr, ringsums[idx])
    print wpeaks

    # fig, ax = plt.subplots()
    # for idx, peak in enumerate(peaks):
    #     plt.plot(lambda_arr, ringsums[idx], label='{0:f}'.format(m1 - idx))
    # ax.legend(title='Order #')
    # ax.set_title("Ring Sums")
    # ax.set_xlabel("Wavelength", fontsize=20)
    # ax.set_ylabel("Counts", fontsize=20)
    # plt.tight_layout()
    # fig.savefig("Argon_with_calibrated_d_L.pdf")
    # plt.show()


    xpadL = np.arange(-n_dlambda / 2., 0) * delta_lambda + lambda_min
    xpadR = np.arange(1, n_dlambda / 2.) * delta_lambda + lambda_max + delta_lambda
    xpad = np.hstack((xpadL, lambda_arr, xpadR))

    instr_width = .007264
    instr = fitting.norm_lorentzian(xpad, 1.0, w, instr_width)
    # instr = instr / instr.sum()  # normalize the kernel

    Ti = 0.6 # eV
    # gamma = 7.7e-5 * w * np.sqrt(Ti/40.0)
    gamma = 3.265e-5 * w * np.sqrt(Ti/40.0)

    g = fitting.gaussian(lambda_arr, 1.0, w, gamma)

    cv = np.convolve(instr, g, mode='valid')

    fig, ax = plt.subplots()
    ax.plot(lambda_arr, ringsums[0], 'r')
    ax1 = plt.twinx()
    ax1.plot(lambda_arr, cv, 'b')
    plt.show()

