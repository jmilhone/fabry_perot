from __future__ import print_function, division
import numpy as np

def bin_data(data, npts=100):
    n, m = divmod(len(data), npts)
    if m != 0:
        d = np.pad(data, pad_width=(0, npts-m), mode='constant', constant_values=np.nan)
        d = np.reshape(d, (n+1, npts))
    else:
        d = data.copy()
        d = np.reshape(d, (n, npts))
    return np.nanmean(d, axis=1), np.nanstd(d, axis=1)
