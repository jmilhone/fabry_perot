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

analyzer = pymultinest.Analyzer(n_params=4, outputfiles_basename="saves/Ar_solver_run4/fp_full_")
stats = analyzer.get_mode_stats()
mode = stats['modes'][0]

mode_vals = mode['maximum a posterior']

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

plt.plot(r**2, sig, 'r')
plt.plot(r**2, vals, 'b')
plt.show()
