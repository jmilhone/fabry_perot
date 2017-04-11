import numpy as np
import matplotlib.pyplot as plt
import analysis.datahelpers as dh
from scipy.optimize import curve_fit, fmin


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
    smooth_points = kwargs.get("smooth_points", 20)
    thres_val = thres * np.max(data)
    up, down = find_peaks(x, data, thres_val, smooth_points=smooth_points)
    npeaks = min(len(up), len(down))
    if npeaks == 0:
        return []
    fits = []
    peaks = []
    for i in range(npeaks):
        vals = data[up[i]:down[i]].copy()
        val_max = data.max()
        vals /= val_max
        xx = x[up[i]:down[i]]
        fit = curve_fit(gaussian, xx, vals, p0=[1.0, np.mean(xx), 0.5 * (xx[-1] - xx[0])])[0]
        fit[0] *= val_max
        fits += [fit]
        peaks += [fit[1]]

    # print "Peak locations: ", peaks

    if plotit:
        fig, ax = plt.subplots()
        ax.plot([x[0], x[-1]], [thres_val] * 2, '--k')
        ax.plot(x, data, lw=1)
        ax.plot(x[up], data[up], 'or', ms=4)
        ax.plot(x[down], data[down], 'oc', ms=4)
        for i, fit in enumerate(fits):
            xx = x[up[i]:down[i]]
            plt.plot(xx, gaussian(xx, *fit), '--g')
        plt.show()
    return peaks


def find_peaks(x, data, maxval, smooth_points=20):
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
    d = dh.smooth(data, smooth_points)
    parse = np.zeros_like(data)
    parse[np.where(d > maxval)] = 1.0

    up = np.where(np.logical_and(parse[1::] == 1.0, parse[0:-1] == 0.0))[0]
    down = np.where(np.logical_and(parse[1::] == 0.0, parse[0:-1] == 1.0))[0]

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
    voight = np.convolve(gaussian_kernel, lor, mode='valid')
    voight = voight * amp + offset
    # print len(voight), data[idx0:idx0+2], voight[idx0:idx0+2]
    chisq = np.sum((data[idx0:idx1] - voight[idx0:idx1])**2 / sigma[idx0:idx1]**2)
    if plotit:
        plt.plot(data)
        plt.plot(voight)
        plt.show()
    return chisq


def argon_Ti_chisq():
    pass


# def argon_chisq(a, x, data, instr_kernel, w0, idx0, idx1, plotit=False):
#     amp = np.exp(a[0])
#     Ti = np.exp(a[1])
#     width = 7.7e-5 * w0 * np.sqrt(Ti/)
