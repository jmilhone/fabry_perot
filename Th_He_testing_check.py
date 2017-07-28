from __future__ import division
import numpy as np
import pymultinest
import He_model2 as model2
from os.path import join
import model
import cPickle as pickle
import model
import matplotlib.pyplot as plt

savedir = "saves/Th_He_testing_solver0/"
basename = "fp_"

analyzer = pymultinest.Analyzer(n_params=6, outputfiles_basename=savedir+basename)

with open("Th_He_data.p", "rb") as inputfile:
    data = pickle.load(inputfile)

rr_He = data['rHe']
data_He = data['sHe']
inds_He = data['inds_He']
#data_He = data_He[inds_He]
#rr_He = rr_He[inds_He]
#sd_He = 0.05 * data_He + 100.

rr_Th = data['rTh']
data_Th = data['sTh']
inds_Th = data['inds_Th']
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
A_He = mode_vals[5]

print L * .004
print d 
print F
wTh = 487.873302
wAr = 487.98634

muAr = 39.948
muTh = 232.03806

Ti_Th = 1000.0 * .025 / 300.0
Ti_Ar = 0.20
Ti_He = 0.857 #9174
He_amplitudes = [0.5471, 0.12806, 1.4311, 2.2558] # 9174

model_Th = model.forward4(rr_Th, L, d, F ,
          [Ti_Th, Ti_Ar], [muTh, muAr], [A_Th, A_Ar], [wTh, wAr]) 

model_He = model2.model_output(rr_He, L, d, Ti_He, He_amplitudes, Q=F) * A_He

r0_He = 3500.0
r0_Th = 3500.0

i_Th = np.argmax(model_Th)
i_He = np.argmax(model_He)

#model_Th *= np.exp(-(rr_Th / r0_Th)**2) / np.exp(-(rr_Th[i_Th] / r0_Th)**2)
#model_He *= np.exp(-(rr_He / r0_He)**2) / np.exp(-(rr_He[i_He] / r0_He)**2)

fig, ax = plt.subplots()
ax.plot(rr_He, data_He, 'r')
ax.plot(rr_He, model_He, 'b')
plt.show(block=False)

fig1, ax1 = plt.subplots()
ax1.plot(rr_Th, data_Th, 'r')
ax1.plot(rr_Th, model_Th, 'b')
plt.show()



