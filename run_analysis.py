from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
from scipy.optimize import fmin
import time
from scipy.special import wofz

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

if __name__ == "__main__":
    binsize = 0.1
    # w = 487.98634
    w = 487.8733 # Th
    n_dlambda = 1536
    delta_lambda = 0.0001

    folder = "Images"
    fname = join(folder, "thorium_ar_5_min_1.nef")
    bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")
    # fname = join(folder, "0015676_000.nef")
    # bg_fname = join(folder, "0015676_001.nef")
    # fname = join(folder, "0015677_000.nef")
    # bg_fname = join(folder, "0015677_001.nef")

    data = im.get_image_data(fname, bg_fname, color='b')
    im.quick_plot(data)
    plt.show()
    print data.shape
    # L = 37562.1081426  # pixels
    # d = .879999589322  # mm

    # new L and d calculation
    # L = 37553.6193101
    # d = 0.87997638795

    # even newer L and d!
    L = 37554.6185682
    d = 0.879904383091

    # zero velocity L
    # L = 37541.57047419231

    # or zero velocity d
    # d = 0.880657

    x0, y0 = (3067.7, 2037.6)

    x0, y0 = rs.locate_center(data, x0, y0, binsize=binsize, plotit=False)

    binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    peaks = fitting.peak_and_fit(binarr, sigarr, thres=0.55, plotit=True)
    r1 = peaks[0]
    r2 = peaks[1]
    m1 = 2.0 * d * 1.0e6 * L / np.sqrt(L ** 2 + r1 ** 2) / w
    m2 = 2.0 * d * 1.0e6 * L / np.sqrt(L ** 2 + r2 ** 2) / w
    print m1 - m2
    marr = [2.0 * d * 1.0e6 * L / np.sqrt(L ** 2 + r ** 2) / w for r in peaks]

    print ""
    print marr
    print np.diff(marr)
    print ""

    print 2.0 / np.sqrt( L**2 + peaks[1]**2)
    print 1.0 / np.sqrt(L**2 + peaks[0]**2) + 1.0/np.sqrt(L**2 + peaks[2]**2)
    costheta = L / np.sqrt(L**2 + r1**2)
    # print (2 * d * 1e6 * costheta / np.floor(2*d*1e6/w) - w)/w * 3.e8
    # costheta = L / np.sqrt(L**2 + (r1+1.0)**2)
    # print (2 * d * 1e6 * costheta / np.floor(2*d*1e6/w) - w)/w * 3.e8
    # print m1, m2
    vc = 1 - 2*d*1.e6/w / np.floor(2*d*1.e6/w) * costheta
    print vc, vc*3.e8
    temp = 1./np.sqrt(1+(r1/L)**2) - 1./np.sqrt(1 + (r2/L)**2)
    print temp
    print 2 * d / w
    v = 1 - 2 * d * 1.e6 * temp / w
    # v *= 3.0e8
    print v, v*3.0e8

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
    print [temp - w for temp in wpeaks]
    print [(temp - w) / w * 3.e8 for temp in wpeaks]
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
    i_widths = [0.00682, 0.00664, 0.00697, 0.00750]

    # instr_width = .0068
    # instr_width = .0067
    # w_instr = np.arange(-512.0 / 2, 512.0 / 2, 1) * 0.0001 + 487.8733

    for jj, instr_width in enumerate(i_widths):
        instr = fitting.norm_lorentzian(xpad, 1.0, w, instr_width)
        # instr = fitting.norm_lorentzian(w_instr, 1.0, 487.8733, instr_width)
        # plt.plot(w_instr, instr)
        # plt.show()
        instr = instr / instr.sum()  # normalize the kernel

        offset_guess = np.mean(ringsums[jj][1:100])
        print offset_guess
        p0 = [np.log(np.max(ringsums[jj]) * 0.015), np.log(1.0), np.log(offset_guess)]
        print np.max(ringsums[jj]) * 0.01
        maxloc = np.argmax(ringsums[jj])
        val = 0.01 * ringsums[jj][maxloc]
        try:
            idx0 = np.max(np.where(ringsums[jj][0:maxloc] < val))
        except ValueError, e:
            print e
            idx0 = 0
        try:
            idx1 = np.min(np.where(ringsums[jj][maxloc:] < val)) + maxloc
        except ValueError, e:
            print e
            idx1 = len(ringsums[jj]) - 1
        # idx0 = 0
        # idx1 = len(ringsums[jj]) - 1
        opt = fmin(fitting.argon_Ti_chisq, p0, args=(lambda_arr, ringsums[jj], instr, w, idx0, idx1))
        print np.exp(opt)
        # fitting.argon_Ti_chisq(opt, lambda_arr, ringsums[jj], instr, w, idx0, idx1, plotit=True)
        fig, ax = plt.subplots()
        ax.plot(lambda_arr, ringsums[jj]-offset_guess, 'b', lw=1)
        Ti = 0.8
        hgamma = instr_width / 2.0
        sigma = 3.265e-5 * w * np.sqrt(Ti / 40.0)
        norm = 1./(sigma * np.sqrt(2 * np.pi))
        z = (lambda_arr - w + 1j*hgamma)/(sigma*np.sqrt(2.0))
        vo = np.real(wofz(z)) * norm
        ax1 = ax.twinx()
        ax1.plot(lambda_arr, vo, '--r', lw=1)
        plt.show()
    # Ti = 0.8 # eV
    # gamma = 7.7e-5 * w * np.sqrt(Ti/40.0)
    # gamma = 3.265e-5 * w * np.sqrt(Ti/40.0)

    # g = fitting.gaussian(lambda_arr, 1.0, w, gamma)
    #
    # cv = np.convolve(g, instr, mode='valid')
    # cv = np.convolve(g, instr, mode='same')

    # fig, ax = plt.subplots()
    # ax.plot(lambda_arr, ringsums[0] / ringsums[0].max(), 'r')
    # ax1 = plt.twinx()
    # ax.plot(lambda_arr, cv / cv.max(), 'b')
    # temp = fitting.norm_lorentzian(lambda_arr, 1.0, w, .0067)
    # ax.plot(lambda_arr, temp/temp.max(), 'c')
    # temp = fitting.norm_lorentzian(lambda_arr, 1.0, w, 0.00725)
    # ax.plot(lambda_arr, temp/temp.max(), 'g')
    # plt.show()
