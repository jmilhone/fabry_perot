from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
from scipy.optimize import minimize_scalar, fmin, curve_fit, differential_evolution
import time
import cPickle as pickle
from analysis.datahelpers import smooth
import multinest_solver as mns


np.seterr(invalid='ignore')  # I got sick of invalid values that happening during minimization
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
global_L_limits = (37000, 38000)
global_d_limits = (.87, .89)

Ar_params = {
    'binsize': 0.1,
    'c_lambda': 487.873302,
    'delta_lambda': 0.0001*4,
    'duc': 0.88,
    'n_dlambda': 512 /4 ,
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

def chisq_d_L(a, radius, weights, pks, dw, w0, ndl=512):

    # L = a[0] * (global_L_limits[1] - global_L_limits[0]) + global_L_limits[0]
    # d = a[1] * (global_d_limits[1] - global_d_limits[0]) + global_d_limits[0]
    L = a[0]
    d = np.exp(a[1])
    m0 = np.floor(2 * d * 1.e6 / w0)
    npks = len(pks)



    rj = [L * np.sqrt((2.e6 * d / ((m0-j)*w0))**2 -1.0) for j in range(npks)]
    print m0,L, d
    print rj
    print pks
    chisq = 0.0
    for idx, r in enumerate(rj):
        chisq += (r - pks[idx])**2 / .02**2
    print chisq,"\n"
    # print L, d, m0
    # rings, warr = rs.ringsum(radius, weights, m0, L, d, pks, dw, ndl=ndl)
    # npks = len(rings) - 1
    # lambda_peaks = []
    # for idx in range(npks):
    #     temp = fitting.peak_and_fit(warr, rings[idx], thres=0.55, plotit=False)
    #     if len(temp) > 0.0:
    #         lambda_peaks += [temp[0]]
    #     else:
    #         print "no peaks\n"
    #         return 1.e3
    # chisq = 0.0
    # for peak in lambda_peaks:
    #     chisq += peak**2 / dw**2
    # # print chisq, "\n"
    return chisq


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
        # print temp
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
    # Ly = (Luc - global_L_limits[0]) / (global_L_limits[1] - global_L_limits[0])
    # dy = (duc - global_d_limits[0]) / (global_d_limits[1] - global_d_limits[0])
    # opt = fmin(minimize_shift_for_d_L, x0=[Ly, dy],
    #            args=(radius, weights, m, pks, w, warr, wmin,
    #                  wmax, dw, ndl), ftol=1.e-6)

    # opt = fmin(chisq_d_L, x0=[Ly, dy],
    #            args=(radius, weights, pks, dw, w, ndl), ftol=1.e-6)
    # Lopt = opt[0] * (global_L_limits[1] - global_L_limits[0]) + global_L_limits[0]
    # dopt = opt[1] * (global_d_limits[1] - global_d_limits[0]) + global_d_limits[0]
    t = time.time()
    mybounds = [(145.0/.0039, 155./.0039), (np.log(.88-.1), np.log(.88+.1))]
    res = differential_evolution(chisq_d_L, bounds=mybounds,
                                 args=(radius, weights, pks, dw, w, ndl), tol=1.e-8)
    Lopt, dopt = res.x
    print res
    dopt = np.exp(dopt)
    # print res.success, res.message
    print mybounds
    print time.time()-t
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
    # f_bg=None
    data = im.get_image_data(f, bgname=f_bg, color='b')
    # im.quick_plot(data)
    # plt.show()
    # data = np.load("100M_Ar.npy")
    # data = np.load("Ar_image.npy")
    # data = np.load("Ar_image.npy") + np.load("Th_image.npy")
    #print data.shape
    # x0 = data.shape[1]/2
    # y0 = data.shape[0]/2
    #print x0, y0, "My artificial center"
    times += [time.time()]
    print "Image done reading, {0} seconds".format(times[-1] - times[-2])
    x0, y0 = rs.locate_center(data, x0, y0, binsize=binsize, plotit=False)
    times += [time.time()]
    print "Center found, {0} seconds".format(times[-1] - times[-2])

    binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    new_binarr = np.concatenate(([0.0], binarr))
    n = len(new_binarr)

    # rarr = np.array([0.5 * (new_binarr[i] + new_binarr[i+1]) for i in range(n-1)])
    rarr = binarr
    # plt.plot(rarr, sigarr)
    # plt.show()
    # plt.plot(binarr, sigarr)
    # plt.show()

    # new_peaks = fitting.peak_and_fit2(rarr, sigarr, thres=0.3, plotit=True, smooth_points=10)
    peaks = fitting.peak_and_fit(rarr, sigarr, thres=0.3, plotit=False, smooth_points=10)

    th_peaks = peaks[0::2]
    ar_peaks = peaks[1::2]



    #print "Overwrote Peaks!"
    #peaksAr = [734.47460093036864, 1166.1300952876536, 1477.9307337048338, 1734.9492968399895]
    #peaksTh = [636.93352662090467, 1106.2487027848206, 1431.5667637305719, 1693.6561601893175]
    #print th_peaks
    #print ar_peaks
    # new_ar_peaks = new_peaks[0::2]
    # new_th_peaks = new_peaks[1::2]
    # peaks = fitting.peak_and_fit(binarr, sigarr, thres=0.65, plotit=False)
    # print "new method peaks: ", new_peaks
    print "Ar Peak locations: ", ar_peaks
    print "Th Peak locations: ", th_peaks
    # print "Ar peaks: ", new_ar_peaks
    # print "Th peaks: ", new_th_peaks
    r1 = peaks[0]
    r2 = peaks[1]

    Luc = find_initial_Lguess(peaks[0], peaks[1], c_lambda, duc)
    # Luc = find_L(peaks)
    # print "new luc", Luc
    # m1 = 2.0 * duc * 1.0e6 * Luc / np.sqrt(Luc ** 2 + r1 ** 2) / c_lambda
    # m2 = 2.0 * duc * 1.0e6 * Luc / np.sqrt(Luc ** 2 + r2 ** 2) / c_lambda

    # Recalculate L and d with m1 and m2 being integers
    # L = np.sqrt(((m2 * r2) ** 2 - (m1 * r1) ** 2) / (m1 ** 2 - m2 ** 2))
    # d m1 * c_lambda * np.sqrt(L ** 2 + r1 ** 2) / (2 * L * 1.e6)

    ny, nx = data.shape
    x = np.arange(0, nx, 1)
    y = np.arange(0, ny, 1)
    xx, yy = np.meshgrid(1. * x - x0, 1. * y - y0)
    R = np.sqrt(xx ** 2 + yy ** 2)

    # Flatten data for simplicity
    R = R.flatten()
    data = data.flatten()
    #save_dir = mns.fix_save_directory("test1")
    save_dir="realpeaks1_.5_err"
    analysis = mns.solve_L_d(save_dir, th_peaks, ar_peaks, L_lim=(149.0/.004, 151.0/.004) )
    L, d, Lsd, dsd = mns.L_d_results(analysis, th_peaks, ar_peaks)

    # Now run with finer resolution
    save_dir = save_dir + "_refined"
    analysis = mns.solve_L_d(save_dir, th_peaks, ar_peaks, L_lim=(L-5.0*Lsd, L+5.0*Lsd) )
    L, d, Lsd, dsd = mns.L_d_results(analysis, th_peaks, ar_peaks)
    m0 = 2.0e6 * d / c_lambda
    # i = np.argsort(R)
    # plt.plot(R[i], data[i])
    # plt.show()

    # print "Overriding L and d"
    # L = 150.1 / 0.0039
    # d = 0.881

    # L = 37554.6185682
    # d = 0.879904383091

    # print L, d
    #L = 38584.6420858
    #d = 0.881218270681
    #m0 = np.floor(2*d*1.e6 / c_lambda)

    

    # print "Testing a new method"
    #maxiter = 10
    ## for i in range(maxiter):
    ##     print i
    ##     L = np.sqrt((r2**2*(m0-1)**2 - r1**2*m0**2)/(2 * m0 - 1))
    ##     darr = [0.5 * (m0 - x) * np.sqrt(L**2 + peaks[x]**2) * c_lambda / (1.e6 *L) for x in range(3)]
    ##     print darr
    ##     d = np.mean(darr)
    ##     m0 = np.floor(2 * d * 1.e6 / c_lambda)
    ##     print m0, L, d, "\n"
    ## L, d = find_d_L(37500.0, duc, R, data, m0, peaks, c_lambda, 0.0, 0.0, 0.0, delta_lambda, ndl=512)
    ##print L, d
    ringsums, lambda_arr = rs.ringsum(R, data, m0, L, d, th_peaks, delta_lambda, ndl=n_dlambda)

    w_peaks = []
    for rsum in ringsums:
        # plt.plot(lambda_arr, rsum)
        # pk = fitting.peak_and_fit(lambda_arr, rsum)
        # w_peaks += pk
        i = np.abs(lambda_arr).argmin() # close enough to the peak
        thres = .5 * rsum[i]
        iL = np.max(np.where(rsum[0:i] < thres))
        iR = np.min(np.where(rsum[i:] < thres)) + i
        # print iL, iR, lambda_arr[iL], lambda_arr[iR]
        inds = range(iL, iR+1)
        pk = np.trapz(lambda_arr[inds] * rsum[inds], x=lambda_arr[inds]) / np.trapz(rsum[inds], x=lambda_arr[inds])
        w_peaks += [pk]
    #    # print pk
    ##     print pk
    ## plt.show()
    print w_peaks
    ##print c_lambda - 487.867
    ##print 487.882 - c_lambda


    ##Find fitting region:
    iL = np.abs((lambda_arr + .006)).argmin()
    iR = np.abs((lambda_arr - .006)).argmin()

    ##print iL, iR
    ind = range(iL, iR)

    Ti = 1000.0 * .025 / 300.0
    sigma = 3.265e-5 * c_lambda * np.sqrt(Ti / lampmu)
    ##print 'sigma', sigma
    
    hwhm = []
    for i, pk in enumerate(w_peaks):
        maxval = np.max(ringsums[i][ind])
        xdata = lambda_arr[ind]
        ydata = ringsums[i][ind] / maxval

        opt, cov = curve_fit(lambda x, a, b: fitting.voigt_profile(x, a, b, sigma, pk),
                             xdata, ydata, p0=[0.01, 0.0036])
        print "Peak {0}: HWHM gamma= {1} nm".format(i, opt[1])
        hwhm.append(opt[1]) 

        voigt = fitting.voigt_profile(lambda_arr, opt[0], opt[1], sigma, pk)
        # voigt = fitting.voigt_profile(lambda_arr, .01, .0036, sigma, 0.0)

        fig, ax = plt.subplots()
        # ax.plot(xdata, ydata, 'b')
        ax.plot(lambda_arr, maxval*voigt, 'r')
        # fig, ax = plt.subplots()
        ax.plot(lambda_arr, ringsums[i], 'b')
        plt.show()

    return L, d, hwhm, x0, y0
    ## ax1 = ax.twinx()
    ## gamma = .0036

    #ii = np.argmin(np.abs(lambda_arr))
    # print ii-iL, iR-ii
    # voigt *= ringsums[0][ii]/voigt.max()
    # ax.plot(lambda_arr, voigt, '--r')
    # plt.show()
    # c_lambda_arr = np.arange(-n_dlambda / 2, n_dlambda / 2, 1) * delta_lambda + c_lambda
    # lambda_min = c_lambda_arr.min()
    # lambda_max = c_lambda_arr.max()
    # print lambda_min, lambda_max
    # npeaks = len(peaks)

    # times += [time.time()]
    # L, d = find_d_L(L, d, R, data, m1, peaks, c_lambda, c_lambda_arr,
    #                 lambda_min, lambda_max, delta_lambda, ndl=n_dlambda)
    # # L = 150.0 / .004
    # # d = 0.88
    # times += [time.time()]
    # marr = [2.0 * d * 1.0e6 * L / np.sqrt(L ** 2 + r ** 2) / c_lambda for r in peaks]
    # print "\nOrder Numbers: ", marr
    # print "Order Differences: ", np.diff(marr), "\n"
    # # print "Overwrote L and d"
    # # L =37561.3832261
    # # d =0.879999606156

    # m1 = 2 * d * 1.e6 * L / np.sqrt(L ** 2 + r1 ** 2) / c_lambda
    # print "\nL, d solver took {0:f} seconds".format(times[-1] - times[-2])
    # print "L={0} pixels".format(L)
    # print "d={0} mm".format(d)
    # print "Highest order m={0}\n".format(m1)
    # ringsums = rs.proper_ringsum(R, data, m1, L, d, peaks, lambda_min,
    #                              lambda_max, delta_lambda, ndl=n_dlambda)
    # wpeaks = []
    # for idx in range(npeaks):
    #     wpeaks += fitting.peak_and_fit(c_lambda_arr, ringsums[idx], thres=0.65, plotit=True)
    # # wpeaks = [wp for wp in wpeaks if np.abs(wp-c_lambda) < 0.05]
    # print "Wavelength peaks (nm) for each order:"
    # print wpeaks, "\n"

    # # gammaT = 7.7e-5 * c_lambda * np.sqrt(Tlamp / lampmu)
    # # gammaT = 3.265e-5 * c_lambda * np.sqrt(Tlamp / lampmu)
    # sigmaT = 3.265e-5 * c_lambda * np.sqrt(Tlamp / lampmu)

#     # xpadL = np.arange(-512 / 2., 0) * delta_lambda + lambda_min
#     # xpadR = np.arange(1, 512 / 2.) * delta_lambda + lambda_max + delta_lambda
#     # xpad = np.hstack((xpadL, c_lambda_arr, xpadR))
#     #
#     # g = fitting.norm_gaussian(c_lambda_arr, 1.0, c_lambda, gammaT)
#     # g = g / np.sum(g)

    # # Ar
    # # idx0 = np.abs(c_lambda_arr - 487.867).argmin()
    # # idx1 = np.abs(c_lambda_arr - 487.882).argmin()
    # # Fake Ar
    # idx0 = np.abs(c_lambda_arr - 487.98).argmin()
    # idx1 = np.abs(c_lambda_arr - 488.02).argmin()

#     # He
#     # idx0 = np.abs(c_lambda_arr - 468.615).argmin()
#     # idx1 = np.abs(c_lambda_arr - 468.625).argmin()

    # widths = []
    # for idx in range(npeaks):
    #     srs = smooth(ringsums[idx], 10)
    #     maxval = srs.max()
    #     vals = ringsums[idx] / maxval
    #     # test_vals = fitting.voigt_profile(c_lambda_arr, sigmaT * np.sqrt(2.0 * np.pi), .005, sigmaT, c_lambda)
    #     # print c_lambda_arr
    #     # print test_vals
    #     # fig, ax = plt.subplots()
    #     # ax.plot(c_lambda_arr, test_vals, '--r')
    #     # ax1 = plt.twinx()
    #     # ax1.plot(c_lambda_arr, vals, 'b')
    #     # plt.show()
    #     opt, cov = curve_fit(
    #         lambda x, amp, a: fitting.voigt_profile(x, amp, a, sigmaT, wpeaks[idx]),# + offset,
    #         c_lambda_arr[idx0:idx1], vals[idx0:idx1], p0=[1.0, .005],# np.mean(vals[1:n_dlambda/10])]
    #     )
    #     print opt
    #     widths.append(opt[1])
    #     fig, ax = plt.subplots()
    #     ax.plot(c_lambda_arr, vals, 'b')
    #     #fitvals = opt[2] + fitting.voigt_profile(c_lambda_arr, opt[0], opt[1], sigmaT, wpeaks[idx])
    #     fitvals = fitting.voigt_profile(c_lambda_arr, opt[0], opt[1], sigmaT, wpeaks[idx])
    #     ax.plot(c_lambda_arr, fitvals, '--r')
    #     plt.show()
    #     # print np.mean(vals[1:50])
    #     # p0 = [np.log(0.1), np.log(5.0*gammaT), np.log(np.mean(vals[1:50]))]
    #     # times += [time.time()]
    #     # opt = fmin(fitting.instr_chisq, p0, args=(xpad, vals, g, c_lambda, idx0, idx1), ftol=1e-5)
    #     # times += [time.time()]
    #     #
    #     # print "Instrument function found, {0} seconds".format(times[-1]-times[-2])
    #     # print np.exp(opt)
    #     # print np.exp(opt[0]), np.exp(opt[1]), np.exp(opt[2])
    #     # fitting.instr_chisq(opt, xpad, vals, g, c_lambda, idx0, idx1, plotit=True)
    #     # widths += [np.exp(opt[1])]

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
    #     "d": d,
    #     "HWHM": widths
    # }
    # with open("calibration_data.p", 'wb') as outfile:
    #     pickle.dump(data_to_save, outfile)

    # return L, d, widths


if __name__ == "__main__":
    binsize = 0.1
    folder = "Images"
    fname = join(folder, "thorium_ar_5_min_1.nef")
    bg_fname = join(folder, "thorium_ar_5_min_1_bg.nef")
    # fname = join(folder, "thorium_5_min_he_3.nef")
    # bg_fname = join(folder, "thorium_5_min_he_bg.nef")
    # center_guess = (3068.39, 2031.85)
    center_guess = (3068.56, 2033.17)
    # He guess
    # center_guess = (3040.05627213, 2024.06787634)

    # L, d, width = run_calibration(fname, bg_fname, center_guess, gas='He')
    # L, d, width = run_calibration(fname, bg_fname, center_guess, gas='Ar')
    L, d, hwhm, x0, y0 = run_calibration(fname, bg_fname, center_guess, gas='Ar')
    # print "FWHM: ", [2 * ww for ww in width]

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
