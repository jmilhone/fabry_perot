from __future__ import division
from matplotlib import rcParams
import numpy as np
import rawpy
import image_helpers as im
import ring_sum as rs
from os.path import join
import matplotlib.pyplot as plt
import plottingtools.core as ptools
import json
import pymultinest
import fitting
#from proper_ring_sum_testing import getR

colors = ptools.tableau20_colors()

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'


def get_data_and_ringsum(f, center_guess, f_bg=None, binsize=1.0, useR=False):
    if f[-3:].lower() == "nef":
        data = im.get_image_data(f, bgname=f_bg, color='b')
    elif f[-3:].lower() == "npy":
        data = np.load(f)
    else:
        print "need a valid file"
        return


    x0, y0 = center_guess
    x0, y0 = rs.locate_center(data, x0, y0, binsize=0.1, plotit=True)
    if useR:
        ny, nx = data.shape
        x = np.arange(0, nx, 1)
        y = np.arange(0, ny, 1)

        xx, yy = np.meshgrid(1.*x-x0, 1.*y-y0)
        R = np.sqrt(xx**2 + yy**2)

        binarr = getR(L, d, 478.98634, norders=6)
        sigarr, _ = np.histogram(R, bins=np.concatenate((np.array([0.]), binarr)), weights=data)
    else:
        binarr, sigarr = rs.quick_ringsum(data, x0, y0, binsize=binsize, quadrants=False)
    return binarr, sigarr

shotnum = 9215
folder = "Images"
#cf = "thorium_ar_5_min_1.nef"
#cf_bg = "thorium_ar_5_min_1_bg.nef"
#f = "0015676_000.nef"
#bg = "0015676_001.nef"
cf = "Th_lamp_488_calib_5m.nef"
cf_bg = None
f = "{0:07d}_000.nef".format(shotnum)
bg = "{0:07d}_001.nef".format(shotnum)

cf = join(folder, cf)
#cf_bg = join(folder, cf_bg)
cf_bg = None
f = join(folder, f)
bg = join(folder, bg)

#folder = "FM_saves"
#cf = "model1.npy" 
#cf = join(folder,cf) 
##f = "Ar_V_3.4_model1.npy"
#f= "Ar_V_neg1.9_ti_.38_model1.npy"
#f = join(folder, f)
#cf_bg=None
#bg=None
#center_guess = (3018., 2010.)
center_guess = (3040., 2004.)
#ccenter_guess = (3018., 2010.)
ccenter_guess = (3071.074026, 2029.44632598)
#center_guess = (3067.04176744, 2035.79588755)

#ccenter_guess = (3040.78197222, 2003.29777455) # old Th calib center guess
#center_guess = (3040.78197222, 2003.29777455) # old Th calib center guess
print f
r, s = get_data_and_ringsum(f, center_guess, f_bg=bg, binsize=1.0, useR=False)
s -= s.min()
#LL, RR = fitting.get_peaks(r, s, thres2=0.05, npks=8)
#inds = [] 
#for i in range(4):
#    inds += range(LL[i], RR[i])
#amp0 = np.max(s[LL[0]:RR[0]])
#print amp0
#pk0 = np.trapz( s[LL[0]:RR[0]] * r[LL[0]:RR[0]], x=r[LL[0]:RR[0]]) / np.trapz( s[LL[0]:RR[0]] , x=r[LL[0]:RR[0]]) 
#print pk0
#print amp0 / np.exp(-(pk0 / 3000.0)**2)
#data = {'r': r.tolist(), "sig": s.tolist(), "L": L, "d": d, "F":F, "idx": inds, "A": amp0 }
#data = {'r': r.tolist(), "sig": s.tolist(), "L": L, "d": d, "F":F, "idx": inds, "A": amp0 / np.exp(-(pk0 / 3000.0)**2)}
#plt.plot(r, s)
#plt.show()
#data = {'r': r.tolist(), "sig": s.tolist(), "idx": inds, "A": amp0 / np.exp(-(pk0 / 3000.0)**2)}

#with open("{0:07d}_data2.json".format(shotnum), 'w') as outfile:
#    print outfile
#    json.dump(data, outfile, indent=4)


#analyzer = pymultinest.Analyzer(n_params=7, outputfiles_basename="saves/solver3_run5/fp_full_")
##analyzer = pymultinest.Analyzer(n_params=8, outputfiles_basename="saves/full_solver_run17/fp_full_")
#stats = analyzer.get_mode_stats()
#local_log_ev = [x['local log-evidence'] for x in stats['modes']]
#ix = np.argmax(local_log_ev)
#mode = stats['modes'][ix]
#
#mode_vals = mode['maximum a posterior']
## 8 param
##L = mode_vals[0]
##d = mode_vals[1]
##F = mode_vals[2]
##Ti0 = mode_vals[3]
##Ti1 = mode_vals[4]
##Amp_Th = mode_vals[5]
##Amp_Ar = mode_vals[6]
##rscale = mode_vals[7]
## 7 param
#L = mode_vals[0]
#d = mode_vals[1]
##print "overwrote L and d"
##L = 150.1 / .004
##d = 0.881
#
#F = mode_vals[2]
##Ti0 = mode_vals[3]
#Ti0 = 1000.0 * .025 / 300.0
#Ti1 = mode_vals[3]
#Amp_Th = mode_vals[4]
#Amp_Ar = mode_vals[5]
#rscale = mode_vals[6]
#
rc, sc = get_data_and_ringsum(cf, ccenter_guess, f_bg=cf_bg, binsize=1.0)
##r, s = get_data_and_ringsum(f, center_guess, f_bg=bg, binsize=1.0, useR=True)
#r, s = get_data_and_ringsum(f, center_guess, f_bg=bg, binsize=1.0, useR=False)
#rr, ss = get_data_and_ringsum(f, center_guess, f_bg=bg, binsize=1.0, useR=False)
##print "Subtracting minimum from argon plasma ring sum"
##s -= s.min()
##ss -= ss.min()
##s -= s[-1]
#
#plt.plot(np.diff(rr), label='CA')
#plt.plot(np.diff(r), label="PRS")
#plt.legend()
#plt.show()
#
##LL, RR = fitting.get_peaks(r, s, thres2=0.1, npks=8)
#LL, RR = fitting.get_peaks(r, s, thres2=0.05, npks=8)
#inds = [] 
#for i in range(4):
#    inds += range(LL[i], RR[i])
#amp0 = np.max(s[LL[0]:RR[0]])
#print amp0
#pk0 = np.trapz( s[LL[0]:RR[0]] * r[LL[0]:RR[0]], x=r[LL[0]:RR[0]]) / np.trapz( s[LL[0]:RR[0]] , x=r[LL[0]:RR[0]]) 
#print pk0
#print amp0 / np.exp(-(pk0 / 3000.0)**2)
##data = {'r': r.tolist(), "sig": s.tolist(), "L": L, "d": d, "F":F, "idx": inds, "A": amp0 }
#data = {'r': r.tolist(), "sig": s.tolist(), "L": L, "d": d, "F":F, "idx": inds, "A": amp0 / np.exp(-(pk0 / 3000.0)**2)}
#with open("{0:07d}_data.json".format(shotnum), 'w') as outfile:
#    print outfile
#    json.dump(data, outfile, indent=4)
#
fig, ax = plt.subplots()
#ax.plot(r**2, s, color=colors[0], label='proper')
##ax1 = ax.twinx()
##ax1.plot(rr**2, ss, color=colors[6])
#ax.plot(rr**2, ss, color=colors[6], label='CA')
##`plt.plot(r[inds]**2, s[inds], color=colors[6])
#plt.legend()
#plt.show()
#fig, ax = plt.subplots()
ax.plot(rc**2, sc, label="Th Lamp", color=colors[0])
ax1 = ax.twinx()
ax1.plot(r**2, s, label="Ar plasma", color=colors[6])
ptools.combine_legend(ax, ax1, frameon=False)
ax.set_xlabel(r"R${}^2$ (px${}^2$)", fontsize=20)
ax.set_ylabel("Counts", fontsize=20, color=colors[0])
ax1.set_ylabel("Counts", fontsize=20, color=colors[6])
ptools.add_thick_box(ax, minor=False)
ptools.add_thick_box(ax1, minor=False)
ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
ax1.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
plt.tight_layout()
plt.show()
