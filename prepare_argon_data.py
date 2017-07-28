from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import ring_sum as rs
import fitting
import json
import image_helpers as im
from os.path import join
from scipy.stats import norm




folder = "synthetic_data/test1/"
#fname = join(folder, "Ar_noV.npy")
#fname = join(folder, "Ar_noV_Ti_1.1.npy")
#fname = join(folder, "Ar_noV_Ti_1.1_smallnoise.npy")
fname = join(folder, "Ar_V_4.7.npy")
#fname = join(folder, "Ar_V_4.7_Ti_1.1.npy")
#fname = join(folder, "Ar_V_4.7_Ti_1.1_smallnoise.npy")
#fname = join(folder, "Ar_V_n3.2_Ti_1.1_smallnoise.npy")
data = np.load(fname)

fig, ax = plt.subplots()
ax.imshow(data, cmap='gray')
ax.set_axis_off()
plt.show()
#bg_data = np.abs(norm(scale=0.05).rvs(data.shape))
bg_data = np.abs(norm(scale=0.2).rvs(data.shape))
x0 = 3018.0
y0 = 2010.0

x0, y0 = rs.locate_center(data, x0, y0, binsize=0.1, plotit=False)
r, sig = rs.quick_ringsum(data, x0, y0, binsize=1.0, quadrants=False)
#sig -= sig.min()
_, bgsig = rs.quick_ringsum(bg_data, x0, y0, binsize=1.0, quadrants=False)
sig -= bgsig

r = np.concatenate(([0.0], r))
r = 0.5*(r[0:-1]+r[1:])
LL, RR = fitting.get_peaks(r, sig, thres2=0.015, npks=4)
inds = []
for L, R in zip(LL, RR):
    i = range(L, R+1)
    inds.extend(i)

outdata = {'r': r.tolist(), "sig": sig.tolist(), "idx": inds, "A": np.max(sig)}
#with open(join(folder,"Ar_noV_data.json"),'w') as outfile:
#with open(join(folder,"Ar_V_4.7_bgfix_data.json"),'w') as outfile:
#with open(join(folder, "Ar_V_4.7_Ti_1.1_data.json"), 'w') as outfile:
#with open(join(folder, "Ar_V_4.7_Ti_1.1_smallnoise_data.json"), 'w') as outfile:
#with open(join(folder, "Ar_V_n3.2_Ti_1.1_smallnoise_bgfix_data.json"), 'w') as outfile:
#with open(join(folder, "Ar_noV_Ti_1.1_data.json"), 'w') as outfile:
#with open(join(folder, "Ar_noV_Ti_1.1_smallnoise_data.json"), 'w') as outfile:
#with open(join(folder, "Ar_noV_Ti_1.1_smallnoise_data_bgfix.json"), 'w') as outfile:
#    json.dump(outdata, outfile, indent=4)

plt.plot(r, sig)
#plt.plot(r, bgsig,'c')
plt.plot(r[inds], sig[inds])
plt.show()
