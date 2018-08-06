from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, RealData
import pymultinest
from ..tools.plotting import my_hist


def determine_fit_range(r, signal, pkr, thres=0.15, plotit=False):
    """Helper function to pick fitting range using peak guess and a threshold

    Args:
        r (np.ndarray): r data
        signal (np.ndarray): ringsum data
        pkr (float): guess for peak location
        thres (float): threshold from peak that is included, default=0.15
        plotit (bools): plots the determined range for checking, default=False

    Returns:
        list: list of indexes within threshold of peak guess
    """
    pkix = np.abs(pkr - r).argmin()
    sig = signal - signal.min()
    pkthr = (1. - thres) * sig[pkix]
    diff = 1.
    i = 1

    Lix = None
    Rix = None

    while diff >= 0.:
        Lix = pkix - i
        diff = sig[Lix] - pkthr
        i += 1
    diff = 1.
    i = 1
    while diff >= 0.:
        Rix = pkix + i
        diff = sig[Rix] - pkthr
        i += 1

    if plotit:
        f, ax = plt.subplots(figsize=(10, 6))
        ax.plot(r, sig, 'b')
        ax.plot(r[Lix - 1:Rix], sig[Lix - 1:Rix], 'r')
        plt.show()

    return range(Lix - 1, Rix + 1)


def find_maximum(x, y, returnval=False):
    """Finds maximum location of a single peak using simple derivative method. Assumes only one
        peak is in given data range.

    Args:
        x (np.ndarray): x array of values
        y (np.ndarray): y array of values
        returnval (bool): Returns the value at the maximum if True, default=False

    Returns:
        pk (float): x position of peak
            if a minimum is detected, returns None
    """
    dy = np.diff(y)
    xx = 0.5 * (x[0:-1] + x[1:])
    i = np.where(np.logical_and(dy[0:-1] > 0.0, dy[1:] < 0.0))[0]
    if len(i) > 1:
        vals = y[i]
        idx = np.argmax(vals)
        i = i[idx]
    val = y[i]
    slope = (dy[i + 1] - dy[i]) / (xx[i + 1] - xx[i])

    if slope < 0.0:
        pk = xx[i] - dy[i] / slope
        if type(pk) is np.ndarray:
            pk_out = pk[0]
        else:
            pk_out = pk
        if returnval:
            return pk_out, val
        else:
            return pk_out
    else:
        return None


def _gauss(beta, x):
    """Gaussian fitting function modified to fit in x^2 space

    y = A*Exp(-0.5*(x^2-x0^2)^2/sigma^4)

    Args:
        beta (Union[list, np.ndarray]): Gaussian parameters
        x (np.ndarray): x points

    Returns:
        np.ndarray: gaussian evaluated at x with parameters beta
    """
    return beta[0] * np.exp(-0.5 * (x ** 2 - beta[1] ** 2) ** 2 / beta[2] ** 4)


def _gauss_fjacb(beta, x):
    """
    derivative of _gauss function wrt to parameters
    dy/dA = y/A
    dy/dx0 = (2*x0*(x^2-x0^2))/sigma^4 * y
    dy/dsigma = (2*(x^2-x0^2)^2)/sigma^5 * y
    """
    _ret = np.concatenate((_gauss(beta, x) / beta[0],
                           (2 * beta[1] * (x ** 2 - beta[1] ** 2) * _gauss(beta, x)) / beta[2] ** 4,
                           (2 * (x ** 2 - beta[1] ** 2) ** 2 * _gauss(beta, x)) / beta[2] ** 5))
    _ret.shape = (3,) + x.shape
    return _ret


def _gauss_fjacd(beta, x):
    """
    derivative of _gauss function wrt x
    dy/dx = (2*x*(x0^2-x^2))/sigma^4 * y
    """
    return (2. * x * (beta[1] ** 2 - x ** 2) * _gauss(beta, x)) / beta[2] ** 4


def find_peak(x, y, x_sd, y_sd, returnval=False, plotit=False):
    """Uses a gaussian fit to find maximum via scipy.odr fit is performed in x^2 space, but all conversions
    are taken care of inside this function

    Args:
        x (np.ndarray): x values for fit
        y (np.ndarray): y values for fit
        x_sd (np.ndarray): standard deviation in x values (typically binsize)
        y_sd (np.ndarray): standard deviation in y values
        returnval (bool): flag for returning fit value at peak location, default=False
        plotit (bool): flag for plotting fit, default=False

    Returns:
        Tuple (float, float, float): peak location, peak location error, peak value
    """

    #y = A*Exp(-0.5*(x^2-x0^2)^2/sigma^4)
    # beta0 = [y.max(),x.mean(),5*(x.max()-x.min())]
    beta0 = [y.max(), np.sqrt(np.mean(x ** 2)), 5.0 * (x.max() - x.min())]
    print(beta0[1])
    if plotit:
        f, ax = plt.subplots()
        ax.errorbar(x, y, xerr=x_sd, yerr=y_sd, fmt='.', color='C0')
        ax.plot(x, _gauss(beta0, x), color='C1')
        ax.set_title('initial guess')
        plt.show()

    data = RealData(x, y, x_sd, y_sd)
    model = Model(_gauss, fjacb=_gauss_fjacb, fjacd=_gauss_fjacd)

    odr = ODR(data, model, beta0)
    output = odr.run()

    pk = output.beta[1]
    pk_err = output.sd_beta[1]
    val = _gauss(output.beta, pk)

    print('Peak {0:.4f} +/- {1:.4f}'.format(pk, pk_err))
    if plotit:
        output.pprint()
        f, ax = plt.subplots()
        ax.errorbar(x, y, xerr=x_sd, yerr=y_sd, fmt='.', color='C0')
        # ax.fill_between(output.xplus,output.y+output.eps,output.y-output.eps,alpha=0.3,color='C1')
        # nvals = 100
        # npts = len(x)
        # yvals = np.zeros((nvals, npts))
        # beta0 = np.random.normal(loc=output.beta[0], scale=output.sd_beta[0], size=nvals)
        # beta1 = np.random.normal(loc=output.beta[1], scale=output.sd_beta[1], size=nvals)
        # beta2 = np.random.normal(loc=output.beta[2], scale=output.sd_beta[2], size=nvals)

        # for i in xrange(nvals):
        #     yvals[i, :] = _gauss([beta0[i], beta1[i], beta2[i]], x)

        # ymean = np.nanmean(yvals, axis=0)
        # ystd = np.nanstd(yvals, axis=0)

        # ax.fill_between(x, output.y+output.eps,output.y-output.eps,alpha=0.3,color='C1')
        err = np.abs(output.eps)
        ax.fill_between(x, _gauss(output.beta, x) + err, _gauss(output.beta, x) - err, alpha=0.6, color='C1')
        # output.pprint()
        # ax.errorbar(x, _gauss(output.beta, x), yerr=output.eps, color='C2')
        # ax.fill_between(x, ymean-ystd, ymean+ystd, alpha=0.3, color='C1')
        ax.axvspan(pk - pk_err, pk + pk_err, color='C2', alpha=0.3)
        plt.show()

    if returnval:
        return pk, pk_err, val
    else:
        return pk, pk_err


def find_maximum_3(x, y, y_sd, basename, plotit=False):
    """
    guass fit using multinest
    """

    A_lim = [0.8 * y.max(), 1.2 * y.max()]
    X0_lim = [x[0], x[-1]]
    S_lim = [0.5 * (x[-1] ** 2 - x[0] ** 2), 3 * (x[-1] ** 2 - x[0] ** 2)]

    def Prior(cube, ndim, nparams):
        cube[0] = cube[0] * (A_lim[1] - A_lim[0]) + A_lim[0]
        cube[1] = cube[1] * (X0_lim[1] - X0_lim[0]) + X0_lim[0]
        cube[2] = cube[2] * (S_lim[1] - S_lim[0]) + S_lim[0]

    def LogLikelihood(cube, ndim, nparams):
        fit = _gauss(cube, x)
        chisq = np.sum((fit - y) ** 2 / y_sd ** 2)
        return -chisq / 2.0

    nparams = 3
    pymultinest.run(LogLikelihood, Prior, nparams,
                    resume=False, outputfiles_basename=basename)

    result = pymultinest.analyse.Analyzer(nparams, outputfiles_basename=basename)

    samples = result.get_equal_weighted_posterior()

    pks = samples[:, 1]
    pk = np.mean(pks)
    pk_err = np.std(pks)

    if plotit:
        f, ax = plt.subplots(1, 2)
        ax[0].errorbar(x, y, yerr=y_sd, fmt='.', color='C0')
        ax[0].axvspan(pk - pk_err, pk + pk_err, color='C1', alpha=0.3)
        my_hist(ax[1], pks)
        plt.show()

    return pk, pk_err


def gaussian_fit(x, y):
    imax = np.argmax(y)
    ymax = y[imax]
    yy = y / ymax

    if any(val < 0.5 for val in y):
        # try finding the fwhm
        left = np.abs(yy[0:imax] - 0.5).argmin()
        right = np.abs(yy[imax:] - 0.5).argmin() + imax
        fwhm = x[right] - x[left]
    else:
        fwhm = np.max(x) - np.min(x)

    guess = [x[imax], fwhm / 2 / np.sqrt(2*np.log(2)) ]
    #fig, ax = plt.subplots()
    #ax.plot(x, yy, 'o')
    #xx = np.linspace(x.min(), x.max(), 100)
    #ax.plot(xx, _gaussian(xx, *guess))
    #plt.show()

    popt, pcov = curve_fit(_gaussian, x, yy, p0=guess, sigma=0.01 * yy, absolute_sigma=True)

    #fig, ax = plt.subplots()
    #ax.plot(x, yy, 'o')
    #ax.plot(xx, _gaussian(xx, *popt))
    #plt.show()

    # print('popt')
    # print(popt)
    # print(np.sqrt(popt[0]))
    # print(pcov)
    pk = np.sqrt(popt[0])
    pk_sd = 0.5 * np.sqrt(pcov[0, 0]) / pk
    return pk, pk_sd


def _gaussian(x, b, c):
    return np.exp(-0.5 * (x - b) ** 2 / c ** 2)
