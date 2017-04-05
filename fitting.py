import numpy as np
import matplotlib.pyplot as plt
import analysis.datahelpers as dh
from scipy.optimize import curve_fit


def gaussian(x, amp, shift, width):
    return amp * np.exp(-(x-shift)**2 / (2.*width)**2)


def peak_and_fit(x, data, thres=0.55, plotit=False, **kwargs):
#     print "Locating Peaks:"
    smooth_points = kwargs.get("smooth_points", 20)
    thres_val = thres*np.max(data)
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
        fit = curve_fit(gaussian, xx, vals, p0=[1.0, np.mean(xx), 0.5*(xx[-1]-xx[0])])[0]
        fit[0] *= val_max
        fits += [fit]
        peaks += [fit[1]]

    # print "Peak locations: ", peaks

    if plotit:
        fig, ax = plt.subplots()
        ax.plot([x[0], x[-1]], [thres_val]*2, '--k')
        ax.plot(x, data, lw=1)
        ax.plot(x[up], data[up], 'or', ms=4)
        ax.plot(x[down], data[down], 'oc', ms=4)
        for i, fit in enumerate(fits):
            xx = x[up[i]:down[i]]
            plt.plot(xx, gaussian(xx, *fit), '--g')
        plt.show()
    return peaks

def find_peaks(x, data, maxval, smooth_points=20):
    d = dh.smooth(data, smooth_points)
    parse = np.zeros_like(data)
    parse[np.where(d > maxval)] = 1.0

    up = np.where(np.logical_and(parse[1::] == 1.0, parse[0:-1] == 0.0))[0]
    down = np.where(np.logical_and(parse[1::] == 0.0, parse[0:-1] == 1.0))[0]

    # assert len(up) == len(down), "Uh oh...Double check the peak finding"
    # if len(up) != len(down):
    #     plt.plot(x, data)
    #     for i in up:
    #         plt.plot(x[i], data[i], 'o')
    #     for i in down:
    #         plt.plot(x[i], data[i], 'o')
    #     plt.show()
    return up, down


