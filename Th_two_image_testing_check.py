from __future__ import division
import numpy as np
import pymultinest
import He_model2 as model2
from os.path import join
import model
import cPickle as pickle
import model
import matplotlib.pyplot as plt
import image_helpers as im
import ring_sum as rs
from scipy.stats import norm
import plottingtools.core as ptools
folder = "Images/"
f_Th = join(folder, "Th_lamp_488_calib_5m.nef")
f_ThHe = join(folder, "Th_lamp_468_calib_5m.nef")
bg_ThHe = join(folder, "Th_lamp_468_calib_5m_background.nef")
bg_Th = join(folder, "Th_lamp_488_calib_5m_background.nef")
data_Th = im.get_image_data(f_Th, bgname=None, color='b')
data_ThHe = im.get_image_data(f_ThHe, bgname=None, color='b')
bg_ThHe = im.get_image_data(bg_ThHe, bgname=None, color='b')
bg_Th = im.get_image_data(bg_Th, bgname=None, color='b')
#folder = "synthetic_data/test1/"
#f_Th =join(folder,"th_lamp_ar.npy")
#f_ThHe = join(folder, "th_lamp_he.npy")
#data_Th = np.load(f_Th)
#data_ThHe = np.load(f_ThHe)


#xThHe = 3043.9579
#yThHe = 2001.8909

#xTh = 3040.87138083
#yTh =  2003.19835754
#xTh = 3040.87138083
#yTh = 2003.19835754
#xTh = 3018.0
#yTh = 2010.0
xTh = 3040.87138083
yTh = 2003.19835754
#plt.plot(data_ThHe[0:2004, 3042])
#plt.show()
#xTh, yTh = rs.locate_center(data_Th, 3040.0, 2005.0, binsize=0.1, plotit=False)
#data_Th -= 150.0
#data_ThHe -= 80.0
rTh, sTh = rs.quick_ringsum(data_Th, xTh, yTh, binsize=1.0, quadrants=False)

#xThHe = 3042.71439254
#yThHe = 2004.44354106
#xThHe = 3018.0
#yThHe = 2010.0
xThHe = 3042.71439254
yThHe = 2004.44354106
#xThHe, yThHe = rs.locate_center(data_ThHe, 3040.0, 2005.0, binsize=0.1, plotit=True)
rThHe, sThHe = rs.quick_ringsum(data_ThHe, xThHe, yThHe, binsize=1.0, quadrants=False)
_, bg_sThHe = rs.quick_ringsum(bg_ThHe, xThHe, yThHe, binsize=1.0, quadrants=False)
_, bg_sTh = rs.quick_ringsum(bg_Th, xThHe, yThHe, binsize=1.0, quadrants=False)

#bg_Th = np.abs(norm(scale=0.2).rvs(data_Th.shape))
#bg_ThHe = np.abs(norm(scale=0.2).rvs(data_ThHe.shape))
_, bg_sig_Th = rs.quick_ringsum(bg_Th, xTh, yTh, binsize=1.0, quadrants=False)
_, bg_sig_ThHe = rs.quick_ringsum(bg_ThHe, xThHe, yThHe, binsize=1.0, quadrants=False)

#sTh -= bg_sig_Th
#sThHe -= bg_sig_ThHe

rTh = np.hstack(([0.0], rTh))
rTh = 0.5 * (rTh[0:-1] + rTh[1:])

rThHe = np.hstack(([0.0], rThHe))
rThHe = 0.5 * (rThHe[0:-1] + rThHe[1:])
#sThHe -= bg_sThHe

sThHe -= sThHe.min()
sTh -= sTh.min()


savedir = "saves/ThAr_ThHe_testing_solver20/"
basename = "fp_"

#analyzer = pymultinest.Analyzer(n_params=6, outputfiles_basename=savedir+basename)
analyzer = pymultinest.Analyzer(n_params=9, outputfiles_basename=savedir+basename)
#analyzer = pymultinest.Analyzer(n_params=9, outputfiles_basename=savedir+basename)

#with open("Th_ThHe_data.p", "rb") as inputfile:
#    data = pickle.load(inputfile)

#rr_ThHe = data['rThHe']
#data_ThHe = data['sThHe']
#inds_ThHe = data['inds_ThHe']
#data_He = data_He[inds_He]
#rr_He = rr_He[inds_He]
#sd_He = 0.05 * data_He + 100.

#rr_Th = data['rTh']
#data_Th = data['sTh']
#inds_Th = data['inds_Th']
#data_Th = data_Th[inds_Th]
#rr_Th = rr_Th[inds_Th]
#sd_Th = 0.01 * data_Th + 100.

stats = analyzer.get_mode_stats()
local_log_ev = [x['local log-evidence'] for x in stats['modes']]
ix = np.argmax(local_log_ev)

mode = stats['modes'][ix]
mode_vals = mode['maximum a posterior']

L = mode_vals[0]
d = mode_vals[1]
F = mode_vals[2]
A_Ar = mode_vals[3]
A_Th = mode_vals[4]
A_ThHe = mode_vals[5]
r0_Th = mode_vals[6]
r0_ThHe = mode_vals[7]

print L * .004
print d 
print F
print r0_Th
print r0_ThHe
wTh = 487.873302
wAr = 487.98634

muAr = 39.948
muTh = 232.03806

Ti_Th = 1000.0 * .025 / 300.0
#Ti_Ar = 0.20
Ti_Ar = mode_vals[8]
print Ti_Ar
# synthetic
#Ti_Ar = 1.0
#Ti_He = 0.857 #9174
#He_amplitudes = [0.5471, 0.12806, 1.4311, 2.2558] # 9174

QQ = (2.0 * F / np.pi)**2
Tii = .025 * 1000.0 / 300.0 
sig1 = 3.276569e-5 * np.sqrt(Tii/232.0) * 488.0
sig2 = 3.276569e-5 * np.sqrt(Tii/232.0) * 468.0 
#norm1 = sig1 * np.sqrt(2.0 * np.pi)
#norm2 = sig2 * np.sqrt(2.0 *  np.pi)
norm1 = 3.0
norm2 = 3.0

model_Th = model.forward4(rTh, L, d, F ,
          [Ti_Th, Ti_Ar], [muTh, muAr], [A_Th, A_Ar], [wTh, wAr]) 
#model_Th += 6 * 4.1e6 / (1 + QQ) * norm1
#print 6 * 4.1e6 / (1 + QQ) * norm1
#model_He = model2.model_output(rr_He, L, d, Ti_He, He_amplitudes, Q=F) * A_He

model_ThHe = model.forward4(rThHe, L, d, F ,
          [Ti_Th,], [muTh,], [A_ThHe,], [468.6195,]) 
#model_ThHe += 9 * 3.6e6 / (1 + QQ) * norm2
#print 9 * 3.6e6 / (1 + QQ) * norm2
#model_ThHe += mode_vals[8]
model_Th *= np.exp(-(rTh / r0_Th)**2)
model_ThHe *= np.exp(-(rThHe / r0_ThHe)**2)

#r0_He = 3500.0
#r0_Th = 3500.0
#
#i_Th = np.argmax(model_Th)
#i_He = np.argmax(model_He)

#model_Th *= np.exp(-(rr_Th / r0_Th)**2) / np.exp(-(rr_Th[i_Th] / r0_Th)**2)
#model_He *= np.exp(-(rr_He / r0_He)**2) / np.exp(-(rr_He[i_He] / r0_He)**2)


#st = analyzer.get_stats()
#sigmas = ['3sigma']
#alphas = [1.0]

#r = np.linspace(0.0, rTh.max(), 2500)
#fig, ax = plt.subplots()
#ax.plot(rTh**2, sTh, 'r')
#for alpha, sig in zip(alphas, sigmas):
#    up = [x[sig][1] for x in st['marginals']]
#    down = [x[sig][0] for x in st['marginals']]
#
#
#    Th_up = model.forward4(r, up[0], up[1], up[2] ,
#           [Ti_Th, Ti_Ar], [muTh, muAr], [up[4], up[3]], [wTh, wAr]) 
#
#    Th_down = model.forward4(r, down[0], down[1], down[2] ,
#           [Ti_Th, Ti_Ar], [muTh, muAr], [down[4], down[3]], [wTh, wAr]) 
#
#    Th_up *= np.exp(-(r / up[6])**2)
#    Th_down *= np.exp(-(r / down[6])**2)
#
#    #ax.fill_between(r**2, Th_down, Th_up, facecolor='b', alpha=alpha)
#    ax.plot(r**2, Th_down, 'b')
#    ax.plot(r**2, Th_up, 'b')
#
#ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
#ax.set_xlabel(r"R (px${}^2$)")
#ax.set_ylabel("Counts")
#plt.show()

with open("Th_ThHe_data.p", 'rb') as infile:
    data = pickle.load(infile)
print data.keys()

def calculate_peaks(L, d, w, norders=4):
    m = 2.e6 * d / w
    m0 = np.floor(m)
    return [L * np.sqrt( m**2 / (m0 - j)**2 - 1.0) for j in range(norders)]
    #return [np.sqrt(2.0) * L * np.sqrt(1.0 - (m0-j) / m) for j in range(norders)]

ww = [487.800942, 488.120436, 487.436431, 487.8733020, 487.98634, 487.934975, ]
for w in ww:
    print "\n",w
    print [x**2 for x  in calculate_peaks(L, d, w)]


www = np.linspace(487, 489, 1000)
chisq = np.zeros_like(www) 
rsq_vals = [169252, 935684]
for idx, w in enumerate(www):
    pks = calculate_peaks(L, d, w, norders=2)
    chisq[idx] = np.sum([(pks[j] - rsq_vals[j])**2 / 100**2 for j in range(2)])
i = np.argmin(chisq)
print www[i]
fig, ax = plt.subplots()
#ax.plot(rr_ThHe[inds_ThHe], data_ThHe[inds_ThHe], 'r')
ax.plot(rThHe**2, sThHe, label='Data')
#ax.plot(rThHe**2, sThHe+46000.0, label='Data')
#ax.plot(rThHe**2, sThHe+150000.0, label='Data')
ax.plot(rThHe**2, model_ThHe, '--',label='Fit')
ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax.set_xlabel(r"R${}^2$ (px${}^2$)", fontsize=20)
ax.set_ylabel("Counts", fontsize=20)
ax.legend(frameon=False, fontsize=18, loc='upper right')
ptools.add_thick_box(ax, minor=False)
fig.tight_layout()
#fname = join(folder, "Plots/", "synth_test1_ThHe_fit.pdf")
fname = join("Plots/", "calib_ThHe_fit.pdf")
#fname = join("Plots/", "calib_ThHe_nofit.pdf")
fig.savefig(fname)
plt.show(block=False)

fig1, ax1 = plt.subplots()
#ax1.plot(rTh**2, sTh, label='Data')
#ax1.plot(rTh**2, sTh+440000, label='Data')
#ax1.plot(rTh**2, sTh+600000, label='Data')
ax1.plot(rTh**2, sTh, label='Data')
ax1.plot(rTh**2, model_Th, '--',label='Fit')
ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ax1.set_xlabel(r"R${}^2$ (px${}^2$)", fontsize=20)
ax1.set_ylabel("Counts", fontsize=20)
ax1.legend(frameon=False, fontsize=18, loc='upper right')
ptools.add_thick_box(ax1, minor=False)
fig1.tight_layout()
fname = join(folder, "Plots/", "synth_test1_Th_fit.pdf")
fname = join("Plots/", "calib_ThAr_fit.pdf")
#fname = join("Plots/", "calib_ThAr_nofit.pdf")
fig1.savefig(fname)
plt.show()



