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

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

analyzer = pymultinest.Analyzer(n_params=4, outputfiles_basename="saves/Ar_solver_run8/fp_full_")
stats = analyzer.get_mode_stats()
mode = stats['modes'][0]

#mode_vals = mode['maximum a posterior']
mode_vals = mode['mean']
mode_sigma = mode['sigma']
print mode['mean']
print mode['sigma']

Ti = mode_vals[0]
V = mode_vals[1]
r0 = mode_vals[2]
A = mode_vals[3]

with open("0015676_data.json", 'r') as infile:
    data = json.load(infile, parse_float=np.float64)

L = data["L"]
d = data["d"]
F = data["F"]
r = np.array(data['r'])
sig = np.array(data['sig'])

mu = 40.0 
w0 = 487.98634
vals = model.forward(r, L, d, F, Ti, mu, w0, nlambda=512, V=V)
vals *= A * np.exp(-(r/r0)**2)

fig, ax = plt.subplots()
ax.plot(r**2, sig)
#plt.errorbar(r**2, sig, yerr=.03*sig+100, color='r')
ax.plot(r**2, vals, '--')
ptools.add_thick_box(ax, minor=False)
ax.set_xlabel(r"R${}^2$ (px${}^2$)", fontsize=20)
ax.set_ylabel("Counts", fontsize=20)
ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
ax.set_title("Ti = {0:4.3f} +/- {1:4.3f} eV".format(Ti, mode_sigma[0]))
plt.tight_layout()
plt.show()
