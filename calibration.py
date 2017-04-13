from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
from scipy.optimize import minimize_scalar, fmin, curve_fit
import time
import cPickle as pickle
from analysis.datahelpers import smooth

np.seterr(invalid='ignore')  # I got sick of invalid values that happening during minimization
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
global_L_limits = (37000, 38000)
global_d_limits = (.87, .89)

Ar_params = {
    'binsize': 0.1,
    'c_lambda': 487.8733,
    'delta_lambda': 0.0001,
    'duc': 0.88,
    'n_dlambda': 512,
    'Tlamp': 1000.0 * .025 / 300.0,  # .025 eV per 300 K
    'lampmu': 232.0,
}
He_params = {
    'binsize': 0.1*4.0,
    'c_lambda': 468.6195,
    'delta_lambda': 0.0001*4.0,
    'duc': 0.88,
    'n_dlambda': 512/4,
    'Tlamp': 1000.0 * .025 / 300.0,  # .025 eV per 300 K
    'lampmu': 232.0,
}


def minimize_shift_for_d_L(a, radius, weights, m0, pks, w,
                           warr, wmin, wmax, dw, ndl=512):
    """
    Calculates chi squared based on residual of wavelength peaks compared to
    the central wavelength w.  It is varying d (the etalon spacing) and L (the
    camera focal length) to get all of the peaks (npks of them) to line up as
    close as they can to the central wavelength w.

    Args:
        a (list): first element is transformed L, second is transformed d
        radius (np.array): flattened array matrix from camera
        weights (np.array): flattened pixel value matrix from camera
        m (float): highest order number (not an integer value)
        npks (int): number of peaks to use starting at order m
        w (float): central wavelength (in nm)
        warr (np.array): wavelength array centered on w
        wmin (float): minimum wavelength in array
        wmax (float): maximum wavelength in array
        dw (float): delta wavelength spacing

    Returns:
        chi squared based on residual of wavelength peaks compared to central wavelength w
    """
    L = a[0] * (global_L_limits[1] - global_L_limits[0]) + global_L_limits[0]
    d = a[1] * (global_d_limits[1] - global_d_limits[0]) + global_d_limits[0]
    m = [2.0 * d * 1.0e6 * L / np.sqrt(L ** 2 + pk ** 2) / w for pk in pks]
    if np.abs(m[0] - m0) > 1:
        print "M got away from us"
        m = 1. * m0
        ringsum = rs.proper_ringsum(radius, weights, m, L, d, pks[0:-1], wmin, wmax, dw, ndl=ndl)
    else:
        ringsum = rs.proper_ringsum(radius, weights, m[0:-1], L, d, pks[0:-1], wmin, wmax, dw, ndl=ndl)
    lambda_peaks = []
    npks = len(pks[0:-1])
    for idx in range(npks):
        temp = fitting.peak_and_fit(warr, ringsum[idx], thres=0.65)
        if len(temp) > 0.0:
            lambda_peaks += [temp[0]]
        else:
            # print "No peaks found"
            return 1.e3
    chisq = 0.0
    for peak in lambda_peaks:
        chisq += (w - peak) ** 2 / dw ** 2
    return chisq


def find_d_L(Luc, duc, radius, weights, m, pks, w, warr, wmin, wmax, dw, ndl=512):
    """
    Calculated L and d from the optimization process for lining up wavelength peaks
    with the central wavelength w

    Args:
        Luc (float): uncalibrated camera focal length (in pixels)
        duc (float): uncalibrated etalon spacing
        radius (np.array): flattened array matrix from camera
        weights (np.array): flattened pixel value matrix from camera
        m (float): highest order number (not an integer value)
        pks (list):  list of peaks to use starting at order m
        w (float): central wavelength (in nm)
        warr (np.array): wavelength array centered on w
        wmin (np.array): minimum wavelength in array
        wmax (np.array): maximum wavelength in array
        dw (float): delta wavelength spacing

    Returns:
        L (float): camera focal length
        d (float): etalon spacing
    """
    Ly = (Luc - global_L_limits[0]) / (global_L_limits[1] - global_L_limits[0])
    dy = (duc - global_d_limits[0]) / (global_d_limits[1] - global_d_limits[0])
    opt = fmin(minimize_shift_for_d_L, x0=[Ly, dy],
               args=(radius, weights, m, pks, w, warr, wmin,
                     wmax, dw, ndl), ftol=1.e-6)

    Lopt = opt[0] * (global_L_limits[1] - global_L_limits[0]) + global_L_limits[0]
    dopt = opt[1] * (global_d_limits[1] - global_d_limits[0]) + global_d_limits[0]

    return Lopt, dopt


def find_initial_Lguess(radius1, radius2, wavelength, d):
    """
    Finds an initial guess for the camera focal length from the first 2 orders

    Args:
        radius1 (float): location of first peak in pixels
        radius2 (float): location of second peak in pixels
        wavelength (float): central wavelength in nm
        d (float): etalon spacing in mm

    Returns:
        L (float): estimate of the camera focal length
    """

    def f(x):
        return np.abs(2. * d * 1.e6 * x / wavelength *
                      (1. / np.sqrt(x ** 2 + radius2 ** 2) -
                       1. / np.sqrt(x ** 2 + radius1 ** 2)) + 1.0)

    opt = minimize_scalar(f, bounds=(1e4, 1e6), method='bounded')

    return opt.x


# def find_L(radii):
#     def f(x):
#         LHS = 2.0 / np.sqrt(x ** 2 + radii[1] ** 2)
#         RHS = 1.0 / np.sqrt(x ** 2 + radii[0] ** 2) + 1.0 / np.sqrt(x ** 2 + radii[2] ** 2)
#         return np.abs(LHS - RHS)

    # opt = minimize_scalar(f, bounds=(149.5 / .004, 150.5 / .004), method='bounded')
    # Lvals = np.linspace(149.5 / .004, 150.5 / .004)
    # fig, ax = plt.subplots(3)
    # ax[0].plot(Lvals, 2.0 / np.sqrt(Lvals ** 2 + radii[1] ** 2) * 1e8)
    # ax[1].plot(Lvals, (1.0 / np.sqrt(Lvals ** 2 + radii[0] ** 2) + 1.0 / np.sqrt(Lvals ** 2 + radii[2] ** 2)) * 1.e8)
    # ax[2].plot(Lvals, f(Lvals) * 1.e8)
    # plt.show()
    # return opt.x


def run_calibration(f, f_bg, center_guess, gas='Ar'):
    times = [time.time()]
    x0, y0 = center_guess

    if gas == 'Ar':
        binsize = Ar_params['binsize']
        c_lambda = Ar_params['c_lambda']
        delta_lambda = Ar_params['delta_lambda']
        duc = Ar_params['duc']
        n_dlambda = Ar_params['n_dlambda']
        Tlamp = Ar_params['Tlamp']
        lampmu = Ar_params['lampmu']
    elif gas == 'He':
        binsize = He_params['binsize']
        c_lambda = He_params['c_lambda']
        delta_lambda = He_params['delta_lambda']
        duc = He_params['duc']
        n_dlambda = He_params['n_dlambda']
        Tlamp = He_params['Tlamp']
        lampmu = He_params['lampmu']
    else:
        print "Did not recognize gas.  Exiting..."
        return None

    data = im.get_image_data(f, f_bg, color='b')
    # im.quick_plot(data)
    # plt.show()

    times += [time.time()]
    print "Image done reading, {0} seconds".format(times[-1] - times[-2])
    x0, y0 = rs.locate_center(data, x0, y0, binsize=binsize, plotit=False)
    times += [time.time()]
    print "Center found, {0} seconds".format(times[-1] - times[-2])

    binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    # plt.plot(binarr, sigarr)
    # plt.show()

    # peaks = fitting.peak_and_fit(binarr, sigarr, thres=0.55, plotit=True)
    peaks = fitting.peak_and_fit(binarr, sigarr, thres=0.65, plotit=False)
    print "Peak locations: ", peaks
    r1 = peaks[0]
    r2 = peaks[1]

    Luc = find_initial_Lguess(peaks[0], peaks[1], c_lambda, duc)
    # Luc = find_L(peaks)
    # print "new luc", Luc
    m1 = 2.0 * duc * 1.0e6 * Luc / np.sqrt(Luc ** 2 + r1 ** 2) / c_lambda
    m2 = 2.0 * duc * 1.0e6 * Luc / np.sqrt(Luc ** 2 + r2 ** 2) / c_lambda

    # Recalculate L and d with m1 and m2 being integers
    L = np.sqrt(((m2 * r2) ** 2 - (m1 * r1) ** 2) / (m1 ** 2 - m2 ** 2))
    d = m1 * c_lambda * np.sqrt(L ** 2 + r1 ** 2) / (2 * L * 1.e6)

    ny, nx = data.shape
    x = np.arange(1, nx + 1, 1)
    y = np.arange(1, ny + 1, 1)
    xx, yy = np.meshgrid(1. * x - x0, 1. * y - y0)
    R = np.sqrt(xx ** 2 + yy ** 2)

    # Flatten data for simplicity
    R = R.flatten()
    data = data.flatten()

    c_lambda_arr = np.arange(-n_dlambda / 2, n_dlambda / 2, 1) * delta_lambda + c_lambda
    lambda_min = c_lambda_arr.min()
    lambda_max = c_lambda_arr.max()
    npeaks = len(peaks)

    times += [time.time()]
    L, d = find_d_L(L, d, R, data, m1, peaks, c_lambda, c_lambda_arr,
                    lambda_min, lambda_max, delta_lambda, ndl=n_dlambda)
    times += [time.time()]
    marr = [2.0 * d * 1.0e6 * L / np.sqrt(L ** 2 + r ** 2) / c_lambda for r in peaks]
    print "\nOrder Numbers: ", marr
    print "Order Differences: ", np.diff(marr), "\n"
    # print "Overwrote L and d"
    # L =37561.3832261
    # d =0.879999606156

    m1 = 2 * d * 1.e6 * L / np.sqrt(L ** 2 + r1 ** 2) / c_lambda
    print "\nL, d solver took {0:f} seconds".format(times[-1] - times[-2])
    print "L={0} pixels".format(L)
    print "d={0} mm".format(d)
    print "Highest order m={0}\n".format(m1)
    ringsums = rs.proper_ringsum(R, data, m1, L, d, peaks, lambda_min,
                                 lambda_max, delta_lambda, ndl=n_dlambda)
    wpeaks = []
    for idx in range(npeaks):
        wpeaks += fitting.peak_and_fit(c_lambda_arr, ringsums[idx], thres=0.65, plotit=True)
    wpeaks = [wp for wp in wpeaks if np.abs(wp-c_lambda) < 0.05]
    print "Wavelength peaks (nm) for each order:"
    print wpeaks, "\n"

    # gammaT = 7.7e-5 * c_lambda * np.sqrt(Tlamp / lampmu)
    # gammaT = 3.265e-5 * c_lambda * np.sqrt(Tlamp / lampmu)
    sigmaT = 3.265e-5 * c_lambda * np.sqrt(Tlamp / lampmu)

    # xpadL = np.arange(-512 / 2., 0) * delta_lambda + lambda_min
    # xpadR = np.arange(1, 512 / 2.) * delta_lambda + lambda_max + delta_lambda
    # xpad = np.hstack((xpadL, c_lambda_arr, xpadR))
    #
    # g = fitting.norm_gaussian(c_lambda_arr, 1.0, c_lambda, gammaT)
    # g = g / np.sum(g)

    # Ar
    # idx0 = np.abs(c_lambda_arr - 487.867).argmin()
    # idx1 = np.abs(c_lambda_arr - 487.882).argmin()

    # He
    idx0 = np.abs(c_lambda_arr - 468.615).argmin()
    idx1 = np.abs(c_lambda_arr - 468.625).argmin()

    widths = []
    for idx in range(npeaks):
        srs = smooth(ringsums[idx], 10)
        maxval = srs.max()
        vals = ringsums[idx] / maxval
        # test_vals = fitting.voigt_profile(c_lambda_arr, sigmaT * np.sqrt(2.0 * np.pi), .005, sigmaT, c_lambda)
        # print c_lambda_arr
        # print test_vals
        # fig, ax = plt.subplots()
        # ax.plot(c_lambda_arr, test_vals, '--r')
        # ax1 = plt.twinx()
        # ax1.plot(c_lambda_arr, vals, 'b')
        # plt.show()
        opt, cov = curve_fit(
            lambda x, amp, a, offset: fitting.voigt_profile(x, amp, a, sigmaT, wpeaks[idx]) + offset,
            c_lambda_arr[idx0:idx1], vals[idx0:idx1], p0=[1.0, .005, np.mean(vals[1:n_dlambda/10])]
        )
        print opt
        widths.append(opt[1])
        fig, ax = plt.subplots()
        ax.plot(c_lambda_arr, vals, 'b')
        fitvals = opt[2] + fitting.voigt_profile(c_lambda_arr, opt[0], opt[1], sigmaT, wpeaks[idx])
        ax.plot(c_lambda_arr, fitvals, '--r')
        plt.show()
        # print np.mean(vals[1:50])
        # p0 = [np.log(0.1), np.log(5.0*gammaT), np.log(np.mean(vals[1:50]))]
        # times += [time.time()]
        # opt = fmin(fitting.instr_chisq, p0, args=(xpad, vals, g, c_lambda, idx0, idx1), ftol=1e-5)
        # times += [time.time()]
        #
        # print "Instrument function found, {0} seconds".format(times[-1]-times[-2])
        # print np.exp(opt)
        # print np.exp(opt[0]), np.exp(opt[1]), np.exp(opt[2])
        # fitting.instr_chisq(opt, xpad, vals, g, c_lambda, idx0, idx1, plotit=True)
        # widths += [np.exp(opt[1])]

    data_to_save = {
        "ringsum": ringsums,
        "peaks": wpeaks,
        "m": m1,
        "delta_lambda": delta_lambda,
        "n_dlambda": n_dlambda,
        "c_lambda_arr": c_lambda_arr,
        "c_lambda": c_lambda,
        "lambda_min": lambda_min,
        "lambda_max": lambda_max,
        "L": L,
        "d": d,
        "HWHM": widths
    }
    with open("calibration_data.p", 'wb') as outfile:
        pickle.dump(data_to_save, outfile)

    return L, d, widths


if __name__ == "__main__":
    binsize = 0.1
    folder = "Images"
    # fname = join(folder, "thorium_ar_5_min_1.nef")
    # bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")
    fname = join(folder, "thorium_5_min_he_3.nef")
    bg_fname = join(folder, "thorium_5_min_he_bg.nef")
    # center_guess = (3068.39, 2031.85)
    center_guess = (3068.56, 2033.17)
    # He guess
    center_guess = (3040.05627213, 2024.06787634)

    L, d, width = run_calibration(fname, bg_fname, center_guess, gas='He')
    print "FWHM: ", [2 * ww for ww in width]

    # data = im.get_image_data(fname, bg_fname, color='b')

    # # Very close guess for center
    # x0, y0 = (3069., 2032.)
    # # Not very close guess
    # # x0, y0 = (3067., 2035.)
    # # REALLY close
    # # x0, y0 = (3069.68678585, 2032.85668627)
    # x0, y0 = rs.locate_center(data, x0, y0, binsize=binsize, plotit=True)
    # #
    # # x0, y0 = (3068.39, 2031.85)  # Cooper's center

    # binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    # peaks = fitting.peak_and_fit(binarr, sigarr, thres=0.55, plotit=False)

    # duc = 0.88
    # c_lambda = 487.8733
    # r1 = peaks[0]
    # r2 = peaks[1]
    # Luc = find_initial_Lguess(peaks[0], peaks[1], c_lambda, duc)
    # m1 = 2.0 * duc * 1.0e6 * Luc / np.sqrt(Luc ** 2 + r1 ** 2) / c_lambda
    # m2 = 2.0 * duc * 1.0e6 * Luc / np.sqrt(Luc ** 2 + r2 ** 2) / c_lambda

    # # Recalculate L and d with m1 and m2 being integers
    # L = np.sqrt(((m2 * r2) ** 2 - (m1 * r1) ** 2) / (m1 ** 2 - m2 ** 2))
    # d = m1 * c_lambda * np.sqrt(L ** 2 + r1 ** 2) / (2 * L * 1.e6)

    # ny, nx = data.shape
    # x = np.arange(1, nx + 1, 1)
    # y = np.arange(1, ny + 1, 1)
    # xx, yy = np.meshgrid(1. * x - x0, 1. * y - y0)
    # R = np.sqrt(xx ** 2 + yy ** 2)

    # # Flatten data for simplicity
    # R = R.flatten()
    # data = data.flatten()

    # # Set up lambda variables
    # delta_lambda = 0.0001  # in nm
    # n_dlambda = 512  # This is Cooper's choice
    # c_lambda_arr = np.arange(-n_dlambda / 2, n_dlambda / 2, 1) * delta_lambda + c_lambda
    # lambda_min = c_lambda_arr.min()
    # lambda_max = c_lambda_arr.max()
    # npeaks = len(peaks)

    # t0 = time.time()
    # # print "Overwrote L and d right before optimization"
    # # L =37561.1867137
    # # d = 0.879999617648
    # L, d = find_d_L(L, d, R, data, m1, npeaks, c_lambda, c_lambda_arr,
    #                 lambda_min, lambda_max, delta_lambda)

    # m1 = 2 * d * 1.e6 * L / np.sqrt(L ** 2 + r1 ** 2) / c_lambda
    # print "\nL, d solver took {0:f} seconds".format(time.time() - t0)
    # print "L={0} pixels".format(L)
    # print "d={0} mm".format(d)
    # print "Highest order m={0}\n".format(m1)

    # ringsums = rs.proper_ringsum(R, data, m1, L, d, peaks, lambda_min,
    #                              lambda_max, delta_lambda)
    # wpeaks = []
    # for idx in range(npeaks):
    #     wpeaks += fitting.peak_and_fit(c_lambda_arr, ringsums[idx])

    # print "Wavelength peaks (nm) for each order:"
    # print wpeaks, "\n"

#     # fig, ax = plt.subplots()
#     # for idx, peak in enumerate(peaks):
#     #     plt.plot(c_lambda_arr, ringsums[idx], label='{0:f}'.format(m1 - idx))
#     # ax.legend(title='Order #')
#     # ax.set_title("Ring Sums")
#     # ax.set_xlabel("Wavelength", fontsize=20)
#     # ax.set_ylabel("Counts", fontsize=20)
#     # plt.show()

# data_to_save = {
#     "ringsum": ringsums,
#     "peaks": wpeaks,
#     "m": m1,
#     "delta_lambda": delta_lambda,
#     "n_dlambda": n_dlambda,
#     "c_lambda_arr": c_lambda_arr,
#     "c_lambda": c_lambda,
#     "lambda_min": lambda_min,
#     "lambda_max": lambda_max,
#     "L": L,
#     "d": d
# }
# with open("calibration_data.p", 'wb') as outfile:
#     pickle.dump(data_to_save, outfile)

# Tlamp = 1000.0 * .025 / 300.0  # .025 eV per kelvin
# lampmu = 232.0
# gammaT = 7.7e-5 * c_lambda * np.sqrt(Tlamp / lampmu)

# xpadL = np.arange(-512 / 2., 0) * delta_lambda + lambda_min
# xpadR = np.arange(1, 512 / 2.) * delta_lambda + lambda_max + delta_lambda
# xpad = np.hstack((xpadL, c_lambda_arr, xpadR))

# g = fitting.norm_gaussian(c_lambda_arr, 1.0, c_lambda, gammaT)
# # plt.plot(c_lambda_arr, g)
# # plt.show()
# g = g / np.sum(g)

# idx0 = np.abs(c_lambda_arr - 487.867).argmin()
# idx1 = np.abs(c_lambda_arr - 487.882).argmin()
# print idx0, idx1
# p0 = [np.log(0.1), np.log(5.0*gammaT), 0.0]
# srs = smooth(ringsums[0], 10)
# maxval = srs.max()
# # plt.plot(c_lambda_arr, ringsums[0]/maxval)
# # plt.show()
# opt = fmin(fitting.instr_chisq, p0, args=(xpad, ringsums[0]/maxval, g, c_lambda, idx0, idx1))

# print np.exp(opt[0]), np.exp(opt[1]), opt[2], gammaT*9
# #
# fitting.instr_chisq(opt, xpad, ringsums[0]/maxval, g, c_lambda, idx0, idx1, plotit=True)

# Check calibration against Argon line
# c_lambda = 487.98634
# delta_lambda = 0.0001  # in nm
# n_dlambda = 512  # This is Cooper's choice
# c_lambda_arr = np.arange(-n_dlambda / 2, n_dlambda / 2, 1) * delta_lambda + c_lambda
# lambda_min = c_lambda_arr.min()
# lambda_max = c_lambda_arr.max()
# npeaks = len(peaks)
#
# allpeaks = fitting.peak_and_fit(binarr, sigarr, thres=0.2, plotit=True)
# r1 = allpeaks[0]
# m1 = 2 * d * 1.e6 * L / np.sqrt(L ** 2 + r1 ** 2) / c_lambda
#
# ringsums = rs.proper_ringsum(R, data, m1, L, d, peaks, lambda_min,
#                              lambda_max, delta_lambda)
# fig, ax = plt.subplots()
# for idx, peak in enumerate(peaks):
#     plt.plot(c_lambda_arr, ringsums[idx], label='{0:f}'.format(m1 - idx))
# ax.legend(title='Order #')
# ax.set_title("Ring Sums")
# ax.set_xlabel("Wavelength", fontsize=20)
# ax.set_ylabel("Counts", fontsize=20)
# plt.show()
