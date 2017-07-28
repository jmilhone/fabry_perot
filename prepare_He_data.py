from __future__ import division
from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
import json
binsize = 1.0
fname = "0009174_000.nef"
bgname = "0009174_001.nef"

folder = "Images/"
fname = join(folder, fname)
bgname = join(folder, bgname)
data = im.get_image_data(fname, bgname=None, color='b')
x0 = 3043.95786886
y0 = 2001.8908683
x0, y0 = rs.locate_center(data, x0, y0, binsize=0.1, plotit=False)
r, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)

r = np.concatenate(([0.0], r))
r = 0.5 * (r[0:-1] + r[1:])

if bgname:
    bgdata = im.get_image_data(bgname, bgname=None, color='b')
    _, bg_sig = rs.quick_ringsum(bgdata, x0, y0, binsize=binsize, quadrants=False)
    sigarr -= bg_sig

fig, ax = plt.subplots()
ax.imshow(data-bgdata, cmap='gray')
ax.set_axis_off()
plt.show()
sigarr -= sigarr.min()

left = [664, 1099, 1383, ]#1635]
right = [927, 1285, 1546,]# 1763]

fig, ax = plt.subplots()
ax.plot(r, sigarr)
idx = []
for L, R in zip(left, right):
    LL = np.abs(r - L).argmin()
    RR = np.abs(r - R).argmin() + 1
    i = range(LL, RR)
    idx.extend(i)
    plt.plot(r[i], sigarr[i])

outdata = {'r': r.tolist(), "sig": sigarr.tolist(), "idx": idx, "A": np.max(sigarr)}
with open("He_0009174_bin1.0+_3order_data.json", 'w') as outfile:
    json.dump(outdata, outfile, indent=4)

plt.show()
