from __future__ import division
from matplotlib import rcParams
import numpy as np
import He_model2 as model2
import json
from os.path import join
import matplotlib.pyplot as plt
import pymultinest
import plottingtools.core as ptools

L = 148.715306908/ .004
d = 0.872201921125
#F = 18.553482898
# from old calibration Ar image only
#F = 22.9272771929
F = 20.0
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

#datafname = "He_0009174_data.json"
datafname = "He_0009174_bin1.0+_1order_data.json"
#datafname = "He_0009174_bin1.0+_3order_data.json"

#savedir = join("saves/","He_solver_run10/") 
savedir = join("saves/","He_solver_run9/") 
basename = "fp_"




with open(datafname,'r') as datafile:
    data = json.load(datafile, parse_float=np.float64)


    rr = np.array(data['r'])
    ss = np.array(data['sig']) 
maxval = ss.max()
#analyzer = pymultinest.Analyzer(n_params=8, outputfiles_basename=savedir+basename)
analyzer = pymultinest.Analyzer(n_params=7, outputfiles_basename=savedir+basename)
stats = analyzer.get_mode_stats()

local_log_ev = [x['local log-evidence'] for x in stats['modes']]
ix = np.argmax(local_log_ev)
mode = stats['modes'][ix]

mode_vals = mode['mean']
mode_sigma = mode['sigma']

Ti = mode_vals[0]
V = mode_vals[1]
#r0 = mode_vals[7]
A = mode_vals[2]
#run 0 and 1
#Amp = [1.0] + mode_vals[3:7]
# run 2
#Amp = mode_vals[3:7] + [1.0]
Amp = mode_vals[2:7]
yy  = [0.2 * 24.0/4.0, 2.2*9.0/4.0, 1.35*30.0/4.0, 5.0*33.0/4.0, 1.0]
print Amp
print yy
print "Ti (eV): ",Ti
print "V (km/s): ", V
#print "r0 (px): ",r0
print "A (Counts): ",A
r = np.linspace(0.0, rr.max(), 5000)
#lin_out = A * model2.model_output(r, L, d, Ti, Amp, V=V, F=F)
lin_out = maxval*model2.model_output(r, L, d, Ti, Amp, V=V, F=F)
#print "set V to zero!!!"
#lin_out = maxval*model2.model_output(r, L, d, Ti, Amp, V=0.0, F=F)
#lin_out *= np.exp(-(r/r0)**2)

#pks = []
#for ww in model2.wavelengths:
#    pks.append(model2.peak_locations(L, d, ww*(1.0 - V/2.998e5)))
fig, ax = plt.subplots()
ax.plot(rr**2, ss, label="Data")
ax.plot(r**2, lin_out, '--', label="Fit")
#for pk in pks:
#    ax.axvline(pk**2, color='k')
ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
ptools.add_thick_box(ax, minor=False)
ax.set_xlabel(r"R${}^2$ (px${}^2$)", fontsize=20)
ax.set_ylabel("Counts", fontsize=20)
ax.set_xlim(3e5,9e5)
#ax.legend(frameon=False, loc='upper right', fontsize=18)
plt.tight_layout()
fname = "He_zoom_fit.pdf"
#fname = "He_nofit.pdf"
fname = join("Plots/", fname)
#fig.savefig(fname)
plt.show()
