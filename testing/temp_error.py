from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
sys.path.append('../')
from core.synthetic_data import full_pattern, forward_model
from core.ringsum import smAng_ringsum, locate_center
from tools.file_io import h5_2_dict, dict_2_h5
from tools.plotting import center_plot, my_hist, my_hist2d

N = 1000
binsize = 0.1

L = np.random.normal(205., 2., size=N)
d = np.random.normal(3.88, 0.005, size=N)
F = np.random.normal(9.,1.,size=N)
w0 = 487.98634
mu = 39.984
amp = 1.0
temp = np.random.uniform(0.05, 0.8, size=N)
v = 0.0

#f,axs = plt.subplots(2,2,figsize=(12,10))
#axs = axs.flatten()
#my_hist(axs[0], L)
#axs[0].set_xlabel('L (mm)')
#my_hist(axs[1], d)
#axs[1].set_xlabel('d (mm)')
#my_hist(axs[2], F)
#axs[2].set_xlabel('F')
#my_hist(axs[3], temp)
#axs[3].set_xlabel('temp (eV)')
#plt.tight_layout()
#plt.show()

sigs = np.zeros((N,10021))
t0 = time.time()
for i,t in enumerate(temp):
    print "{0} of {1}".format(i+1,N)
    data = full_pattern(L[i], d[i], F[i], w0, mu, temp[i], 0.0, nprocs=24)
    binarr,sig = smAng_ringsum(data, 3007.5, 2007.5, binsize=binsize, quadrants=False)
    binarr = np.concatenate(([0.0],binarr))
    binarr = 0.5 * (binarr[0:-1] + binarr[1:])
    sigs[i,:] = sig
    dict_2_h5('./Data/synthetic/temp_error_test_01.h5', {'L':L, 'd':d, 'F':F, 'w0':w0, 'mu':mu, 'amp':amp, 'temp':temp, 'v':v, 'sigs':sigs, 'binarr':binarr})
print 'It took {0} seconds to do {1} iterations'.format(time.time()-t0,N)     
