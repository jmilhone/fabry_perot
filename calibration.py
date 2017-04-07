from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
from scipy.optimize import minimize_scalar, fmin
import time
import cPickle as pickle
from analysis.datahelpers import smooth

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
global_L_limits = (37000, 38000)
global_d_limits = (.87, .89)


def minimize_shift_for_d_L(a, radius, weights, m, npks, w,
                           warr, wmin, wmax, dw):
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

    ringsum = rs.proper_ringsum(radius, weights, m, L, d, peaks, wmin, wmax, dw)
    lambda_peaks = []
    for idx in range(npks):
        temp = fitting.peak_and_fit(warr, ringsum[idx], thres=0.65)
        if len(temp) > 0.0:
            lambda_peaks += [temp[0]]
        else:
            print "No peaks found"
            return 1.e3
    # print lambda_peaks
    chisq = 0.0
    for peak in lambda_peaks:
        chisq += (w - peak) ** 2 / .0001 ** 2  # binsize is 0.0001 nm
    # print chisq
    return chisq


def find_d_L(Luc, duc, radius, weights, m, npks, w, warr, wmin, wmax, dw):
    """
    Calculated L and d from the optimization process for lining up wavelength peaks
    with the central wavelength w

    Args:
        Luc (float): uncalibrated camera focal length (in pixels)
        duc (float): uncalibrated etalon spacing
        radius (np.array): flattened array matrix from camera
        weights (np.array): flattened pixel value matrix from camera
        m (float): highest order number (not an integer value)
        npks (int): number of peaks to use starting at order m
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
               args=(radius, weights, m, npks, w, warr, wmin,
                     wmax, dw))

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


if __name__ == "__main__":
    binsize = 0.1
    folder = "Images"
    fname = join(folder, "thorium_ar_5_min_1.nef")
    bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")
    data = im.get_image_data(fname, bg_fname, color='b')

    # Very close guess for center
    # x0, y0 = (3069., 2032.)
    # Not very close guess
    # x0, y0 = (3067., 2035.)
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

    # Set up lambda variables
    delta_lambda = 0.0001  # in nm
    n_dlambda = 512  # This is Cooper's choice
    c_lambda_arr = np.arange(-n_dlambda / 2, n_dlambda / 2, 1) * delta_lambda + c_lambda
    lambda_min = c_lambda_arr.min()
    lambda_max = c_lambda_arr.max()
    npeaks = len(peaks)

    t0 = time.time()
    print "Overwrote L and d right before optimization"
    L =37561.1867137
    d = 0.879999617648
    # L, d = find_d_L(L, d, R, data, m1, npeaks, c_lambda, c_lambda_arr,
    #                 lambda_min, lambda_max, delta_lambda)

    m1 = 2 * d * 1.e6 * L / np.sqrt(L ** 2 + r1 ** 2) / c_lambda
    print "\nL, d solver took {0:f} seconds".format(time.time() - t0)
    print "L={0} pixels".format(L)
    print "d={0} mm".format(d)
    print "Highest order m={0}\n".format(m1)

    ringsums = rs.proper_ringsum(R, data, m1, L, d, peaks, lambda_min,
                                 lambda_max, delta_lambda)
    wpeaks = []
    for idx in range(npeaks):
        wpeaks += fitting.peak_and_fit(c_lambda_arr, ringsums[idx])

    print "Wavelength peaks (nm) for each order:"
    print wpeaks, "\n"

    # fig, ax = plt.subplots()
    # for idx, peak in enumerate(peaks):
    #     plt.plot(c_lambda_arr, ringsums[idx], label='{0:f}'.format(m1 - idx))
    # ax.legend(title='Order #')
    # ax.set_title("Ring Sums")
    # ax.set_xlabel("Wavelength", fontsize=20)
    # ax.set_ylabel("Counts", fontsize=20)
    # plt.show()

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
        "d": d
    }
    with open("calibration_data.p", 'wb') as outfile:
        pickle.dump(data_to_save, outfile)

    Tlamp = 1000.0 * .025 / 300.0  # .025 eV per kelvin
    lampmu = 232.0
    gammaT = 7.7e-5 * c_lambda * np.sqrt(Tlamp / lampmu)

    xpadL = np.arange(-512 / 2., 0) * delta_lambda + lambda_min
    xpadR = np.arange(1, 512 / 2.) * delta_lambda + lambda_max + delta_lambda
    xpad = np.hstack((xpadL, c_lambda_arr, xpadR))

    g = fitting.norm_gaussian(c_lambda_arr, 1.0, c_lambda, gammaT)
    # plt.plot(c_lambda_arr, g)
    # plt.show()
    g = g / np.sum(g)

    idx0 = np.abs(c_lambda_arr - 487.867).argmin()
    idx1 = np.abs(c_lambda_arr - 487.882).argmin()
    print idx0, idx1
    p0 = [np.log(0.1), np.log(5.0*gammaT), 0.0]
    rs = smooth(ringsums[0], 10)
    maxval = rs.max()
    # plt.plot(c_lambda_arr, ringsums[0]/maxval)
    # plt.show()
    opt = fmin(fitting.instr_chisq, p0, args=(xpad, ringsums[0]/maxval, g, c_lambda, idx0, idx1))

    print np.exp(opt[0]), np.exp(opt[1]), opt[2], gammaT*9

    fitting.instr_chisq(opt, xpad, ringsums[0]/maxval, g, c_lambda, idx0, idx1, plotit=True)