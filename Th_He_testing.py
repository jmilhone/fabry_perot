from __future__ import division
from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
import cPickle as pickle


folder = "Images/"
f_Th = join(folder, "Th_lamp_488_calib_5m.nef")

f_He = join(folder, "0009174_000.nef")
f_He_bg = join(folder, "0009174_001.nef")

data_Th = im.get_image_data(f_Th, bgname=None, color='b')
data_He = im.get_image_data(f_He, bgname=None, color='b')
bg_He = im.get_image_data(f_He_bg, bgname=None, color='b')


xHe = 3043.9579
yHe = 2001.8909

xTh = 3040.87138083
yTh =  2003.19835754
#xTh, yTh = rs.locate_center(data_Th, 3040.0, 2005.0, binsize=0.1, plotit=False)
rTh, sTh = rs.quick_ringsum(data_Th, xTh, yTh, binsize=1.0, quadrants=False)

rHe, sHe = rs.quick_ringsum(data_He, xHe, yHe, binsize=1.0, quadrants=False)
_, He_bg = rs.quick_ringsum(bg_He, xHe, yHe, binsize=1.0, quadrants=False)

sHe -= He_bg

sHe -= sHe.min()
sTh -= sTh.min()

rmin_He = 670
rmax_He = 970
imin_He = np.abs(rHe - rmin_He).argmin()
imax_He = np.abs(rHe - rmax_He).argmin()
inds_He = range(imin_He, imax_He+1)

rmin_Th = 620
rmax_Th = 760
imin_Th = np.abs(rTh - rmin_Th).argmin()
imax_Th = np.abs(rTh - rmax_Th).argmin()
inds_Th = range(imin_Th, imax_Th+1)

output_data = {'rHe': rHe, 'sHe': sHe, 'rTh': rTh, 'sTh': sTh, 'inds_He': inds_He, 'inds_Th': inds_Th}
with open("Th_He_data.p", 'wb') as outfile:
    pickle.dump(output_data, outfile)

fig, ax = plt.subplots()
ax.plot(rTh, sTh, 'b')
ax.plot(rTh[inds_Th], sTh[inds_Th], 'g')
ax1 = ax.twinx()
ax1.plot(rHe, sHe, 'r')
ax1.plot(rHe[inds_He], sHe[inds_He], 'c')
plt.show()


