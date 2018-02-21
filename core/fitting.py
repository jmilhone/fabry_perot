import numpy as np
import matplotlib.pyplot as plt

def determine_fit_range(r, sig, pkr, thres=0.15, plotit=False):
    '''
    helper function to pick fitting
    range using peak guess and a threshold

    Args:
        r (np.ndarray): r data
        sig (np.ndarray): ringsum data
        pkr (float): guess for peak location
        thresh (float, default=0.15): threshold
            from peak that is included
        plotit (bools, default=False): plots
            the determined range for checking

    Returns:
        idxs (list): list of indexes within threshold
            of peak guess
    '''
    pkix = np.abs(pkr-r).argmin()
    pkthr = (1. - thres) * sig[pkix]
    diff = 1.
    i = 1
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
        f,ax = plt.subplots(figsize=(10,6))
        ax.plot(r, sig, 'b')
        ax.plot(r[Lix-1:Rix], sig[Lix-1:Rix],'r')
        plt.show()

    return range(Lix-1, Rix) 

def find_maximum(x, y):
    '''
    finds maximum location of a single peak using
    simple derivative method. Assumes only one
    peak is in given data range.

    Args:
        x (np.ndarray): x array of values
        y (np.ndarray): y array of values

    Returns:
        pk (float): x position of peak
            if a minimum is detected, returns None
    '''
    dy = np.diff(y)
    xx = 0.5*(x[0:-1] + x[1:])
    i = np.where(np.logical_and(dy[0:-1]>0.0, dy[1:]<0.0))[0]
    if len(i) > 1:
        vals = y[i]
        idx = np.argmax(vals)
        i = i[idx]
    slope = (dy[i+1] - dy[i]) / (xx[i+1] - xx[i])

    if slope < 0.0:
        pk = xx[i] - dy[i] / slope
        if type(pk) is np.ndarray:
            return pk[0]
        else:
            return pk
    else:
        return None

