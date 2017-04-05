from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
from scipy.optimize import minimize_scalar, fmin
import time

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
global_L_limits = (37000, 38000)
global_d_limits = (.87, .89)


def minimize_shift_for_d_L(a, R, weights, m, npeaks, c_lambda, lambda_arr, lambda_min, lambda_max, delta_lambda):
    L = a[0] * (global_L_limits[1]-global_L_limits[0]) + global_L_limits[0]
    d = a[1] * (global_d_limits[1]-global_d_limits[0]) + global_d_limits[0]
    # print L, d
    ringsum = rs.proper_ringsum(R, weights, m, L, d, peaks, lambda_min, lambda_max, delta_lambda)
    lambda_peaks = []
    for idx in range(npeaks):
        temp = fitting.peak_and_fit(lambda_arr, ringsum[idx], thres=0.65)
        if len(temp) > 0.0:
            lambda_peaks += [temp[0]]
        else:
            print "No peaks found"
            return 1.e3
    # print lambda_peaks
    chisq = 0.0
    for peak in lambda_peaks:
        chisq += (c_lambda - peak)**2 / .0001**2  # binsize is 0.0001 nm
    # print chisq
    return chisq


def find_d_L(Luc, duc, radius, weights, m, npks, w, warr, wmin, wmax, dw):

    Ly = (Luc - global_L_limits[0]) / (global_L_limits[1] - global_L_limits[0])
    dy = (duc - global_d_limits[0]) / (global_d_limits[1] - global_d_limits[0])
    opt = fmin(minimize_shift_for_d_L, x0=[Ly, dy],
               args=(radius, weights, m, npks, w, warr, wmin,
                     wmax, dw))

    Luc = opt[0] * (global_L_limits[1] - global_L_limits[0]) + global_L_limits[0]
    duc = opt[1] * (global_d_limits[1] - global_d_limits[0]) + global_d_limits[0]

    return Luc, duc


def find_initial_Lguess(r1, r2, wavelength, d):
    f = lambda x: np.abs(2.*d*1.e6*x/wavelength * (1./np.sqrt(x**2 + r2**2) - 1./np.sqrt(x**2+r1**2)) + 1.0)
    # (d in mm)*1.e6 = d in nm

    opt = minimize_scalar(f, bounds=(1e4, 1e6), method='bounded')

    # L = np.logspace(4, 6, 10000)
    # plt.plot(L, f(L))
    # plt.xlim(10000, 50000)
    # plt.show()
    # print opt
    # print opt.x
    return opt.x

if __name__ == "__main__":
    binsize = 0.1
    folder = "Images"
    fname = join(folder, "thorium_ar_5_min_1.nef")
    bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")
    data = im.get_image_data(fname, bg_fname, color='b')

    #Very close guess for center
    # x0, y0 = (3069., 2032.)
    # REALLY close
    x0, y0 = (3069.68678585, 2032.85668627)
    x0, y0 = rs.locate_center(data, x0, y0, binsize=binsize, plotit=False)
    #
    # x0, y0 = (3068.39, 2031.85)  # Cooper's center

    binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    peaks = fitting.peak_and_fit(binarr, sigarr, thres=0.55, plotit=False)

    duc = 0.88
    c_lambda = 487.8733
    r1 = peaks[0]
    r2 = peaks[1]
    Luc = find_initial_Lguess(peaks[0], peaks[1], c_lambda, duc)
    m1 = 2.0 * duc * 1.0e6 * Luc / np.sqrt(Luc**2 + r1**2) / c_lambda
    m2 = 2.0 * duc * 1.0e6 * Luc / np.sqrt(Luc**2 + r2**2) / c_lambda

    # Recalculate L and d with m1 and m2 being integers
    L = np.sqrt(((m2*r2)**2 - (m1*r1)**2)/(m1**2 - m2**2))
    d = m1 * c_lambda * np.sqrt(L**2 + r1**2) / (2*L*1.e6)


    ny, nx = data.shape
    x = np.arange(1, nx+1, 1)
    y = np.arange(1, ny+1, 1)
    xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
    R = np.sqrt(xx**2 + yy**2)

    # Flatten data for simplicity
    R = R.flatten()
    data = data.flatten()

    # Set up lambda variables
    delta_lambda = 0.0001  # in nm
    n_dlambda = 512  # This is Cooper's choice
    c_lambda_arr = np.arange(-n_dlambda/2, n_dlambda/2, 1)*delta_lambda + c_lambda
    lambda_min = c_lambda_arr.min()
    lambda_max = c_lambda_arr.max()
    npeaks = len(peaks)


    t0 = time.time()
    L, d = find_d_L(L, d, R, data, m1, npeaks, c_lambda, c_lambda_arr,
                    lambda_min, lambda_max, delta_lambda)

    m1 = 2 * d * 1.e6 * L / np.sqrt(L**2 + r1**2) / c_lambda
    print "\nL, d solver took {0:f} seconds".format(time.time()-t0)
    print "L={0} pixels".format(L)
    print "d={0} mm".format(d)
    print "Highest order m={0}\n".format(m1)

    ringsums = rs.proper_ringsum(R, data, m1, L, d, peaks, lambda_min,
                                 lambda_max, delta_lambda)
    wpeaks = []
    for idx in range(npeaks):
        wpeaks += fitting.peak_and_fit(c_lambda_arr, ringsums[idx])

    print "Wavelength peaks (nm) for each order:"
    print wpeaks,"\n"

    fig, ax = plt.subplots()
    for idx, peak in enumerate(peaks):
        plt.plot(c_lambda_arr, ringsums[idx], label='{0:f}'.format(m1-idx))
    ax.legend(title='Order #')
    ax.set_title("Ring Sums")
    ax.set_xlabel("Wavelength", fontsize=20)
    ax.set_ylabel("Counts", fontsize=20)
    plt.show()

    # ringsums = rs.proper_ringsum(R, data, m1, L, d, peaks, lambda_min, lambda_max, delta_lambda)
    # old_peaks = []
    # for i in range(npeaks):
    #     old_peaks += fitting.peak_and_fit(c_lambda_arr, ringsums[i])
    # # find_d_L([L, d], R, data, m1, peaks, c_lambda, c_lambda_arr, lambda_min, lambda_max,
    # #          delta_lambda)
    # Ly = (L-global_L_limits[0]) / (global_L_limits[1]-global_L_limits[0])
    # dy = (d-global_d_limits[0]) / (global_d_limits[1]-global_d_limits[0])
    # t0 = time.time()
    # opt = fmin(minimize_shift_for_d_L, x0=[Ly, dy],
    #            args=(R, data, m1, peaks, c_lambda, c_lambda_arr, lambda_min,
    #                  lambda_max, delta_lambda))
    # print time.time()-t0

    # fig, ax = plt.subplots()
    # for idx, peak in enumerate(peaks):
    #     plt.plot(c_lambda_arr, ringsums[idx], label='{0:f}'.format(m1))
    # ax.legend()
    # ax.set_title("Inital L, d ringsum")
    # plt.show(block=False)

    # print "Intial", L, d, m1
    # L = opt[0] * (global_L_limits[1]-global_L_limits[0]) + global_L_limits[0]
    # d = opt[1] * (global_d_limits[1]-global_d_limits[0]) + global_d_limits[0]
    # m1 = 2 * d * 1.e6 * L / np.sqrt(L**2 + r1**2) / c_lambda

    # ringsums_final = rs.proper_ringsum(R, data, m1, L, d, peaks, lambda_min, lambda_max, delta_lambda)
    # peaks = []
    # for i in range(npeaks):
    #     peaks += fitting.peak_and_fit(c_lambda_arr, ringsums_final[i])#, thres=.35)
    # print "Final", L, d, m1
    # print old_peaks
    # print peaks
    # fig, ax = plt.subplots()
    # print len(c_lambda_arr)
    # for idx, peak in enumerate(peaks):
    #     plt.plot(c_lambda_arr, ringsums_final[idx], label='{0:f}'.format(m1-idx))
    # ax.legend()
    # ax.set_title("Final L, d ringsum")
    # plt.show()
    # print time.time() - t1
    # for i in range(npeaks):
    #     t2 = time.time()
    #     m = m1 - i
    #     rmin = np.sqrt((2*d*L*1.e6/(m*lambda_max))**2 - L**2)
    #     rmax = np.sqrt((2*d*L*1.e6/(m*lambda_min))**2 - L**2)
    #     inds = np.where(np.logical_and(R > rmin, R<rmax))[0]

        # r0 = 2.*d*1.e6*L / m / delta_lambda
        # # print r0, rmin, rmax

        # dr = np.sqrt((1./(1./np.sqrt(L**2 + rmin**2) - 1./r0))**2 - L**2) - rmin
        # rR = [rmin, rmin + dr]
        # j = 0
        # while rR[-1] < rmax:
        #     j += 1
        #     dr = np.sqrt((1. / (1. / np.sqrt(L ** 2 + rR[-1] ** 2) - 1. / r0)) ** 2 - L ** 2) - rR[-1]
        #     # print j, dr, rR[-1], rmax
        #     rR += [rR[-1] + dr]

#         # print len(rR), len(c_lambda_arr)

        # sig, _ = np.histogram(R[inds], bins=rR, weights=data[inds])
        # # to correspond to my calibration wavelength array, sig is actually backwards...
        # sig = sig[::-1]
        # print time.time() - t2
        # plt.plot(c_lambda_arr, sig)
        # plt.title("idx: {0} order: {1}".format(i, m))
        # plt.show()

