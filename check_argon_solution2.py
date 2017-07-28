from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
#import forward_model as fm
import time
import plottingtools.core as ptools
import cPickle as pickle
import pymultinest
import model
import json
from calibration_solver import read_L_d_results, read_finesse_results
from os.path import join
shotnum = 9215
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

#analyzer = pymultinest.Analyzer(n_params=4, outputfiles_basename="saves/solver_Ar4/fp_full_")
#analyzer = pymultinest.Analyzer(n_params=4, outputfiles_basename="saves/modifiedArsolver_run1/fp_full_")
#analyzer = pymultinest.Analyzer(n_params=4, outputfiles_basename="saves/Ar_solver_syntest1_17/fp_full_")
analyzer = pymultinest.Analyzer(n_params=4, outputfiles_basename="saves/Ar_solver_run25/fp_full_")
stats = analyzer.get_mode_stats()

local_log_ev = [x['local log-evidence'] for x in stats['modes']]
ix = np.argmax(local_log_ev)
mode = stats['modes'][ix]

#mode_vals = mode['maximum a posterior']
mode_vals = mode['mean']
mode_sigma = mode['sigma']
print mode['mean']
print mode['sigma']

Ti = mode_vals[0]
V = mode_vals[1]
r0 = mode_vals[2]
A = mode_vals[3]
print "Ti (eV):", Ti
print "V (km/s):",  V
print "r0 (px):",r0
#with open("0015676_data.json", 'r') as infile:
#with open("{0:07d}_data.json".format(shotnum), 'r') as infile:
with open("{0:07d}_data2.json".format(shotnum), 'r') as infile:
#datafname = join("synthetic_data/test1/", "Ar_V_4.7_data.json")
#datafname = join("synthetic_data/test1/", "Ar_V_4.7_bgfix_data.json")
#datafname = join("synthetic_data/test1/", "Ar_noV_Ti_1.1_data.json")
#datafname = join("synthetic_data/test1/", "Ar_noV_Ti_1.1_smallnoise_data_bgfix.json")
#datafname = join("synthetic_data/test1/", "Ar_V_4.7_Ti_1.1_data.json")
#datafname = join("synthetic_data/test1/", "Ar_V_4.7_Ti_1.1_smallnoise_data.json")
#datafname = join("synthetic_data/test1/", "Ar_V_n3.2_Ti_1.1_smallnoise_data.json")
#datafname = join("synthetic_data/test1/", "Ar_V_n3.2_Ti_1.1_smallnoise_bgfix_data.json")
#datafname = join("synthetic_data/test1/", "Ar_noV_data.json")
#with open(datafname, 'r') as infile:
    data = json.load(infile, parse_float=np.float64)
#L, d = read_L_d_results("saves/Ld_test11")
#F = read_finesse_results("saves/finesse_solver7")


#F = 19.77149120640
#L = 148.599558491 / .004
#d = 0.872203104696
#L =148.715985377 / .004
#d = 0.872201612198
#F = 20.2484267098
#L = data["L"]
#d = data["d"]
#L = 149.973834053 / .004
#d = 0.877040367054
#F = 22.2018789225

# real calib with fake F
L = 148.715306908/ .004
d = 0.872201921125
F = 20.0
#print "overwrote L and d  and v to be perfect model"
#L = 150.1 / .004
#d = .881
#V = -1.8
#F = data["F"]
r = np.array(data['r'])
sig = np.array(data['sig'])
rr = np.linspace(0.0, r.max(), 10000)
#mu = 40.0 
mu = 39.948
w0 = 487.98634
#vals = model.forward(r, L, d, F, Ti, mu, w0, nlambda=512, V=V)
vals = model.forward3(rr, L, d, F, Ti, mu, w0, nlambda=512, V=V)
vals *= A * np.exp(-(rr/r0)**2)
#vals *= A * np.exp(-(r/r0)**2)

fig, ax = plt.subplots()
ax.plot(r**2, sig, label="Data")
idx = data['idx']
#ax.plot(r[idx]**2, sig[idx], 'g')
#plt.errorbar(r**2, sig, yerr=.03*sig+100, color='r')
#ax.plot(r**2, vals, '--')
#ax.plot(rr**2, vals, '--', label="Fit")
ptools.add_thick_box(ax, minor=False)
ax.set_xlabel(r"R${}^2$ (px${}^2$)", fontsize=20)
ax.set_ylabel("Counts", fontsize=20)
ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
#ax.set_title("Ti = {0:4.3f} +/- {1:4.3f} eV".format(Ti, mode_sigma[0]))
#ax.legend(frameon=False, loc='upper right', fontsize=18)
plt.tight_layout()
#fname = "Ar_9215_fit.pdf"
fname = "Ar_9215_nofit.pdf"
fname = join("Plots/",fname)
fig.savefig(fname)
plt.show()
