import numpy as np
import matplotlib.pyplot as plt

def my_hist(ax, data, bins=None):
    if bins is not None:
        hist, bins = np.histogram(data, density=True, bins=bins)
    else:
        hist, bins = np.histogram(data, density=True)

    bw = bins[1]-bins[0]

    ax.bar(bins[0:-1], hist*bw, width=bw)


