from __future__ import division
from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from os.path import join
import fitting
import cPickle as pickle
from scipy.stats import norm

folder = "Images/"
f_Th = join(folder, "Th_lamp_488_calib_5m.nef")
#
f_ThHe = join(folder, "Th_lamp_468_calib_5m.nef")
bg_ThHe = join(folder, "Th_lamp_468_calib_5m_background.nef")
bg_Th = join(folder, "Th_lamp_488_calib_5m_background.nef")
#folder = "synthetic_data/test1/"
#f_Th =join(folder,"th_lamp_ar.npy")
#f_ThHe = join(folder, "th_lamp_he.npy")
#data_Th = np.load(f_Th)
#data_ThHe = np.load(f_ThHe)

data_Th = im.get_image_data(f_Th, bgname=None, color='b')
data_ThHe = im.get_image_data(f_ThHe, bgname=None, color='b')
bg_ThHe = im.get_image_data(bg_ThHe, bgname=None, color='b')
bg_Th = im.get_image_data(bg_Th, bgname=None, color='b')

#fig, ax = plt.subplots()
#ax.imshow(data_ThHe, cmap='gray')
#ax.set_axis_off()
#plt.show()
#xThHe = 3043.9579
#yThHe = 2001.8909

#xTh = 3040.87138083
#yTh =  2003.19835754
#xTh, yTh = 3018.0, 2010.0
#data_Th -= 150.0
#data_ThHe -= 80.0

xTh, yTh = rs.locate_center(data_Th, 3040.0, 2005.0, binsize=0.1, plotit=False)
rTh, sTh = rs.quick_ringsum(data_Th, xTh, yTh, binsize=1.0, quadrants=False)

#xThHe, yThHe = 3018.0, 2010.0
xThHe, yThHe = rs.locate_center(data_ThHe, 3040.0, 2005.0, binsize=0.1, plotit=True)
rThHe, sThHe = rs.quick_ringsum(data_ThHe, xThHe, yThHe, binsize=1.0, quadrants=False)
#_, bg_sThHe = rs.quick_ringsum(bg_ThHe, xThHe, yThHe, binsize=1.0, quadrants=False)

#bg_Th = np.abs(norm(scale=0.2).rvs(data_Th.shape))
#bg_ThHe = np.abs(norm(scale=0.2).rvs(data_ThHe.shape))

_, bg_sig_Th = rs.quick_ringsum(bg_Th, xTh, yTh, binsize=1.0, quadrants=False)
_, bg_sig_ThHe = rs.quick_ringsum(bg_ThHe, xThHe, yThHe, binsize=1.0, quadrants=False)
#plt.plot(rThHe, sThHe)
#plt.plot(rThHe, bg_sig_ThHe)
#plt.show()
#sTh -= bg_sig_Th
#sThHe -= bg_sig_ThHe

rTh = np.hstack(([0.0], rTh))
rTh = 0.5 * (rTh[0:-1] + rTh[1:])

rThHe = np.hstack(([0.0], rThHe))
rThHe = 0.5 * (rThHe[0:-1] + rThHe[1:])

#sThHe -= bg_sThHe

sThHe -= sThHe.min()
sTh -= sTh.min()
# From previous analysis, looking for a 1/(1+Q) offset value
sThHe += 150000.0
sTh += 600000.0
rmin_ThHe = 540.0
rmax_ThHe = 597.0
#rmin_ThHe = 614.0 
#rmax_ThHe = 654.0
imin_ThHe = np.abs(rThHe - rmin_ThHe).argmin()
imax_ThHe = np.abs(rThHe - rmax_ThHe).argmin()
inds_ThHe = range(imin_ThHe, imax_ThHe+1)

new_ThHe = None
rrThHe = None
lefts_ThHe = [540.0, 1018.0, 1334.0, 1587.0]
rights_ThHe = [597.0, 1044.0, 1352.0, 1606.0]
# synthetic_data/test1/
#lefts_ThHe = [156, 873, 1237, 1511]
#rights_ThHe = [326, 923, 1264, 1532]

maxThHe = np.max(sThHe)
maxTh = np.max(sTh)

for L, R in zip(lefts_ThHe, rights_ThHe):
    i1 = np.abs(rThHe - L).argmin()
    i2 = np.abs(rThHe - R).argmin()
    i = range(i1,i2+1)
    temp = sThHe[i] 
    #temp /= temp.max()
    #temp *= maxThHe
    if new_ThHe is None:
        new_ThHe = temp.copy()
        rrThHe = rThHe[i].copy()
    else:
        new_ThHe = np.hstack( (new_ThHe, temp.copy()))
        rrThHe = np.hstack( (rrThHe, rThHe[i].copy()))

new_Th = None
rrTh = None
#lefts_Th = [620.0, 1077.0, 1389.0]
#rights_Th = [760.0, 1159.0, 1456.0]

lefts_Th = [625.0, 1080.0, 1393.0, 1646.0]
rights_Th = [755.0, 1155.0, 1453.0, 1700.0]

#synthetic_data/test1/
#lefts_Th = [502, 1018, 1343, 1610]
#rights_Th = [686, 1121, 1429, 1680]
for L, R in zip(lefts_Th, rights_Th):
    i1 = np.abs(rTh - L).argmin()
    i2 = np.abs(rTh - R).argmin()
    i = range(i1,i2+1)
    temp = sTh[i] 
    #temp /= temp.max()
    #temp *= maxTh
    if new_Th is None:
        new_Th = temp.copy()
        rrTh = rTh[i].copy()
    else:
        new_Th = np.hstack( (new_Th, temp.copy()))
        rrTh = np.hstack( (rrTh, rTh[i].copy()))
print sThHe.min()
plt.plot(rThHe, sThHe, 'r')
plt.plot(rrThHe, new_ThHe, 'b')
plt.show()

plt.plot(rTh, sTh, 'r')
plt.plot(rrTh, new_Th, 'b')
plt.show()

rmin_Th = 620
rmax_Th = 760
imin_Th = np.abs(rTh - rmin_Th).argmin()
imax_Th = np.abs(rTh - rmax_Th).argmin()
inds_Th = range(imin_Th, imax_Th+1)

#output_data = {'rThHe': rThHe, 'sThHe': sThHe, 'rTh': rTh, 'sTh': sTh, 'inds_ThHe': inds_ThHe, 'inds_Th': inds_Th}
#with open("Th_ThHe_data.p", 'wb') as outfile:
#    pickle.dump(output_data, outfile)

output_data = {'rThHe': rrThHe, 'sThHe': new_ThHe, 'rTh': rrTh, 'sTh': new_Th }
#with open(join(folder,"Th_ThHe_data_bgfix.p"), 'wb') as outfile:
with open("Th_ThHe_data_bghack3_realdata.p", 'wb') as outfile:
    pickle.dump(output_data, outfile)

#fig, ax = plt.subplots()
#ax.plot(rTh, sTh, 'b')
#ax.plot(rTh[inds_Th], sTh[inds_Th], 'g')
#ax1 = ax.twinx()
#ax1.plot(rThHe, sThHe, 'r')
#ax1.plot(rThHe[inds_ThHe], sThHe[inds_ThHe], 'c')
#plt.show()


