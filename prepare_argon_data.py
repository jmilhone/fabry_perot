from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import ring_sum as rs
import fitting
import json
import image_helpers as im
from os.path import join
from scipy.stats import norm
import sys
import fp_helpers



shotnum = int(sys.argv[1])

folder = "Images/" 
out_folder = "Calibration/Ti_saves/"
out_folder = join(out_folder, "{0:07d}".format(shotnum))
fp_helpers.make_directory(out_folder)
#folder = "synthetic_data/test1/"
#fname = join(folder, "Ar_noV.npy")
#fname = join(folder, "Ar_noV_Ti_1.1.npy")
#fname = join(folder, "Ar_noV_Ti_1.1_smallnoise.npy")
#fname = join(folder, "Ar_V_4.7.npy")
#fname = join(folder, "Ar_V_4.7_Ti_1.1.npy")
#fname = join(folder, "Ar_V_4.7_Ti_1.1_smallnoise.npy")
#fname = join(folder, "Ar_V_n3.2_Ti_1.1_smallnoise.npy")

fname = join(folder, "{0:07d}_000.nef".format(shotnum))
data = im.get_image_data(fname, color='b')
bgfname = join(folder, "{0:07d}_001.nef".format(shotnum))
bg_data = im.get_image_data(bgfname, color='b')
#data = np.load(fname)

fig, ax = plt.subplots()
ax.imshow(data, cmap='gray')
ax.set_axis_off()
plt.show()
#bg_data = np.abs(norm(scale=0.05).rvs(data.shape))
#bg_data = np.abs(norm(scale=0.2).rvs(data.shape))
x0 =2880.95466832 
y0 = 2028.22974993

x0, y0 = rs.locate_center(data, x0, y0, binsize=0.1, plotit=False)
r, sig = rs.quick_ringsum(data, x0, y0, binsize=0.1, quadrants=False)
#sig -= sig.min()
_, bgsig = rs.quick_ringsum(bg_data, x0, y0, binsize=0.1, quadrants=False)
plt.plot(r, sig, 'b')
plt.plot(r, bgsig, 'g')
plt.show()
sig -= bgsig

r = np.concatenate(([0.0], r))
r = 0.5*(r[0:-1]+r[1:])
LL, RR = fitting.get_peaks(r, sig, thres2=0.025, npks=4)
inds = []
for L, R in zip(LL, RR):
    i = range(L, R+1)
    inds.append(i)

outdata = {'r': r.tolist(), "sig": sig.tolist(), "idx": inds, "A": np.max(sig)}
#with open(join(out_folder, "{0:07d}_data.json".format(shotnum)), 'w') as outfile:
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
i, j = inds[0][0], inds[0][-1]
x = range(i-50, j+50)
plt.plot(r[x], sig[x], 'g')
#for x in inds:
#    plt.plot(r[x], sig[x])
plt.show()
