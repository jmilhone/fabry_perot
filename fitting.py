import numpy as np
import matplotlib.pyplot as plt
import analysis.datahelpers as dh
from scipy.optimize import curve_fit, fmin
from scipy.special import wofz

def gaussian(x, amp, shift, width):
    """
    Returns guassian evaluated on x

    amp * exp( -(x-shit)**2/width**2)

    Args:
        x (np.array): x values to evaluate on
        amp (float): amplitude of gaussian
        shift (float): x shift of guassian
        width (float): shift of gaussian

    Returns:
        np.array or float (depends on x)
    """
    return amp * np.exp(-(x - shift) ** 2 / (2. * width) ** 2)

def peak_and_fit2(x, data, thres=0.55, plotit=False, **kwargs):
    smooth_points = kwargs.get("smooth_points", 5)
    thres_val = thres * np.max(data)
    up, down = find_peaks(x, data, thres_val, smooth_points=smooth_points, plotit=True)

    peaks2fit = kwargs.get('npeaks', 10)

    npeaks = min(len(up), len(down))
    if npeaks == 0:
        return []
    else:
        npeaks = min(npeaks, peaks2fit)
    if plotit:
        fig, ax = plt.subplots()
        ax.plot(x**2, data)
    peaks = []
    for i in range(npeaks):
        vals = data[up[i]:down[i]].copy()
        new_thres = 0.8 * np.max(vals)
        xx = x[up[i]:down[i]].copy()

        pk =  np.sqrt(np.trapz(xx**2 * vals, x=xx**2) / np.trapz(vals, x=xx**2))
        peaks.append(pk)
        # print up[i] - down[i]
        # u, d = find_peaks(xx, vals, new_thres, smooth_points=10)
        # new_fig = plt.figure()
        # plt.plot(xx[u[0]:d[0]]**2, dh.smooth(vals[u[0]:d[0]], 5))
        # plt.show()
        # inds = range(u[0], d[0])
        # new_vals = dh.smooth(vals[inds], 1)
        # new_vals /= new_vals.max()
        # new_xx = xx[inds]**2
        # fit, fit_cov = curve_fit(gaussian, new_xx, new_vals,
        #                          p0=[1.0, np.mean(new_xx), 0.5*(new_xx[-1] - new_xx[0])])

        # print u, d
        # print xx[u[0]], xx[d[0]], np.sqrt(fit[1])

        # if plotit:
        #     plt.plot(xx[u[0]]**2, vals[u[0]], 'or')
        #     plt.plot(xx[d[0]]**2, vals[d[0]], 'og')
    if plotit:
        plt.show()
    return peaks

def peak_and_fit(x, data, thres=0.55, plotit=False, **kwargs):
    """
    Find peaks in the data and returns a list of x locations

    Args:
        x (np.array): x values for data
        data (np.array): y values for x
        thres (float): look for peaks above thres * max(data)
        plotit (bool): plot if True
        **kwargs:
            smooth_points (int): number of points to smooth
    Returns:
        peaks (list): locations of the peaks in x
    """
    smooth_points = kwargs.get("smooth_points", 5)
    thres_val = thres * np.max(data)
    up, down = find_peaks(x, data, thres_val, smooth_points=smooth_points)
    peaks2fit = kwargs.get('npeaks', 10)
    npeaks = min(len(up), len(down))
    if npeaks == 0:
        return []
    else:
        npeaks = min(npeaks, peaks2fit)
    fits = []
    peaks = []
    for i in range(npeaks):
        vals = data[up[i]:down[i]].copy()
        val_max = data.max()
        vals /= val_max
        xx = x[up[i]:down[i]]**2
        # print np.trapz(xx * vals, x=xx) / np.trapz(vals, x=xx)
        # print type(vals), type(xx)
        # print xx
        # print vals
        # plt.plot(xx, vals)
        # plt.show()
        if len(xx) > 5:
            # new_thres = .3 * np.max(data[up[i]:down[i]])
            # L = np.max(np.where(data[0:up[i]] < new_thres))
            # R = np.min(np.where(data[down[i]:] < new_thres)) + down[i]

            fit, fit_cov = curve_fit(gaussian, xx, vals, p0=[1.0, np.mean(xx), 0.5 * (xx[-1] - xx[0])])
            # fit, fit_cov = curve_fit(gaussian, x[L:R], data[L:R] / data[L:R].max(), p0=[1.0, np.mean(x[L:R]), 0.5 * (x[R] - x[L])])
            # print fit_cov
            fit[0] *= val_max
            fits += [fit]
            peaks += [np.sqrt(fit[1])]


            # print "\n",L, up[i]
            # print R, down[i]

            # print np.trapz(x[L:R] * data[L:R], x=x[L:R]) / np.trapz(data[L:R], x=x[L:R])
    # print "Peak locations: ", peaks
    if plotit:
        fig, ax = plt.subplots()
        ax.plot([x[0]**2, x[-1]**2], [thres_val] * 2, '--k')
        ax.plot(x**2, data, lw=1)
        ax.plot(x[up]**2, data[up], 'or', ms=4)
        ax.plot(x[down]**2, data[down], 'oc', ms=4)
        for i, fit in enumerate(fits):
            xx = x[up[i]:down[i]]**2
            plt.plot(xx, gaussian(xx, *fit), '--g')
        ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        ax.set_xlabel(r"$r^2$ (Pixels${}^2$)")
        ax.set_ylabel("Counts")
        plt.tight_layout()
        plt.show()
    # if plotit:
    #     fig, ax = plt.subplots()
    #     ax.plot([x[0], x[-1]], [thres_val] * 2, '--k')
    #     ax.plot(x, data, lw=1)
    #     ax.plot(x[up], data[up], 'or', ms=4)
    #     ax.plot(x[down], data[down], 'oc', ms=4)
    #     for i, fit in enumerate(fits):
    #         xx = x[up[i]:down[i]]
    #         plt.plot(xx, gaussian(xx, *fit), '--g')
    #     plt.show()
    # print peaks
    return peaks


def find_peaks(x, data, maxval, smooth_points=5, plotit=False):
    """

    Args:
        x (np.array): x values for data
        data (np.array): y values for x
        maxval (float): only consider points in data above maxval
        smooth_points (int): number of points used for smoothing

    Returns:
        up (list): list of the x values that start on the upward side of peak
            at maxval
        down (list): list of the x values that end on the downward side of
            peak at maxval
    """
    if smooth_points > 1:
        d = dh.smooth(data, smooth_points)
    else:
        d = data.copy()

    parse = np.zeros_like(data)
    parse[np.where(d >= maxval)] = 1.0

    up = np.where(np.logical_and(parse[1::] == 1.0, parse[0:-1] == 0.0))[0]
    down = np.where(np.logical_and(parse[1::] == 0.0, parse[0:-1] == 1.0))[0]
    if plotit:
        plt.plot(x,d)
        for xx in up:
            plt.plot(x[xx], d[xx], 'oc')
        for xx in down:
            plt.plot(x[xx], d[xx], 'og')
        plt.show()
    return up, down

def norm_gaussian(x, a, shift, gamma):
    """
    Returns a normalized gaussian at shift with amplitude a and width gamma.  The area is a.

    G(x) = a * np.exp(-0.5*((x-shift)/gamma)**2) / gamma / np.sqrt(2 * np.pi)
    Args:
        x (np.array): x values
        a (float): amplitude
        shift (float): center of gaussian
        gamma (float): width of gaussian
    """
    return a * np.exp(-0.5*((x-shift)/gamma)**2) / gamma / np.sqrt(2 * np.pi)


def norm_lorentzian(x, a, shift, gamma):
    """
    Returns a normalized lorentzian at shift with amplitude a and width gamma.  The area is a.

    L(x) = a* 0.5*gamma / ((x-shift)**2 + (0.5*gamma)**2) / np.pi
    Args:
        x (np.array): x values
        a (float): amplitude
        shift (float): center of gaussian
        gamma (float): width of gaussian
    """
    half_gamma = 0.5*gamma
    return a * half_gamma / ((x-shift)**2 + half_gamma**2) / np.pi


def instr_chisq(a, x, data, gaussian_kernel, w0, idx0, idx1, plotit=False):
    """
    Calculates chi squared with the gaussian kernel convolved with
    lorentzian instrument function with the a parameters compared to
    the ringsum data in the index range of (idx0, idx1)

    Args:
        a (list): parameters for finding instrument function
            (amp, width, offset)
        x (np.array): wavelength array to eval instrument function on
        data (np.array): ring sum to match to
        gaussian_kernel (np.array): lamp doppler broadenned gaussian
            Note: sum(gaussian_kernel) == 1 !!!
        w0 (float): wavelength shift in nm
        idx0 (int): left index for fitting region
        idx1 (int): right index for fitting region

    Returns:
        chisq (float)
    """
    amp = np.exp(a[0])
    width = np.exp(a[1])
    offset = np.exp(a[2])
    lor = norm_lorentzian(x, 1.0, w0, width)
    sigma = np.sqrt(data) + 10.0  #10 is to avoid division by zero
    voigt = np.convolve(gaussian_kernel, lor, mode='valid')
    voigt = voigt * amp + offset
    # print len(voigt), data[idx0:idx0+2], voigt[idx0:idx0+2]
    chisq = np.sum((data[idx0:idx1] - voigt[idx0:idx1])**2 / sigma[idx0:idx1]**2)
    if plotit:
        plt.plot(data, lw=1)
        plt.plot(voigt, '--', lw=1)
        plt.show()
    return chisq


def argon_Ti_chisq(a, x, data, lor_kernel, w0, idx0, idx1, plotit=False):
    """

    Args:
        a:
        x:
        data:
        lor_kernel:
        w0:
        idx0:
        idx1:
    """
    amp = np.exp(a[0])
    Ti = np.exp(a[1])
    offset = np.exp(a[2])
    width = 3.265e-5 * w0 * np.sqrt(Ti / 40.0)
    g = norm_gaussian(x, 1.0, w0, width)

    voigt = np.convolve(g, lor_kernel, mode='valid')
    # voigt = np.convolve(g, lor_kernel, mode='same')
    voigt = amp * voigt + offset
    n = len(data[idx0:idx1])
    sigma = np.sqrt(np.abs(data)) + 10.0
    chisq = np.sum((data[idx0:idx1] - voigt[idx0:idx1])**2 / sigma[idx0:idx1]**2) / (n - 4.0)
    if plotit:
        print np.trapz(g, x=x)
        fig, ax = plt.subplots()
        ax.plot(x, data, 'r', lw=1)
        # ax1 = ax.twinx()
        # ax1.plot(voigt, '--b')
        ax.plot(x, voigt, 'b', lw=1)
        low, high = ax.get_ylim()

        left = x[idx0]
        right = x[idx1]

        ax.plot([left]*2, [low, high], '--k')
        ax.plot([right]*2, [low, high], '--k')

        # h = np.max(data) / np.max(lor_kernel)
        # npts = len(x)
        # print npts, len(lor_kernel)
        # ax.plot(x+.03, h*lor_kernel[npts/2:npts+npts/2], '--g', lw=1)
        plt.tight_layout()
        plt.show()
    return chisq


def voigt_profile(x, amp, gamma, sigma, mu):
    """
    G(x;sigma,mu) = exp(-(x-mu)**2/(2 sigma)**2) / (sigma sqrt(2 pi))
    L(x;gamma,mu) = gamma/pi / ( (x-mu)**2 + gamma**2)
    gamma is half-width at half maximum for L(x;gamma,mu)

    V(x;sigma,gamma,mu) = Re[w(z)]/(sigma sqrt(2pi))
    z = (x-mu + i*gamma) / (sigma sqrt(2))

    Returns:

    """
    z = x - mu + 1j*gamma
    z *= 1.0 / (sigma * np.sqrt(2.0))
    norm = sigma * np.sqrt(2.0 * np.pi)
    return amp * np.real(wofz(z)) / norm


# def argon_chisq(a, x, data, instr_kernel, w0, idx0, idx1, plotit=False):
#     amp = np.exp(a[0])
#     Ti = np.exp(a[1])
#     width = 7.7e-5 * w0 * np.sqrt(Ti/)
