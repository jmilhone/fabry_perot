from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import plottingtools.core as ptools
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from os.path import join
import json
from histogram import my_hist
import fp_helpers
import argparse
import model as model
import cPickle as pickle
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

w = 487.98634
mu = 39.948

shotnum = 26248
folder = join("Calibration/Ti_saves", "{0:07d}".format(shotnum))
#folder = join("Calibration/Ti_saves", "{0:07d}_1".format(shotnum))
Ti_post, A_post, V_post = fp_helpers.read_Ar_Ti_results(folder)

Lpost, dpost = fp_helpers.read_Ld_results("Calibration/Ld_saves/2017_09_25/")
Fpost, _, _, _ = fp_helpers.read_finesse_results("Calibration/Finesse_saves/2017_09_25")

data_fname = "{0:07d}_data.json".format(shotnum)
data_fname = join(folder, data_fname)
with open(data_fname, 'r') as datafile:
    data = json.load(datafile, parse_float=np.float64)

order = 0
r = np.array(data['r'])
s = np.array(data['sig'])
idx = data['idx'][order]




fig, ax = plt.subplots()
my_hist(ax, Ti_post, bins=15)
ax.set_xlabel('Ti (eV)')
ax.set_ylabel("P(Ti)")
plt.show(block=False)

fig, ax = plt.subplots()
my_hist(ax, A_post, bins=15)
ax.set_xlabel('A (Counts)')
ax.set_ylabel("P(A)")
plt.show(block=False)

fig, ax = plt.subplots()
my_hist(ax, V_post, bins=15)
ax.set_xlabel('V (km/s)')
ax.set_ylabel("P(V)")
plt.show(block=True)

rr = r[idx]
print len(Ti_post)
Ti_post = Ti_post[::12]
A_post = A_post[::12]
V_post = V_post[::12]

Lpost = Lpost[::150]
dpost = dpost[::150]
Fpost = Fpost[::25]

print Ti_post.shape, Lpost.shape, Fpost.shape

x = len(Lpost)
y = len(Fpost)
z = len(Ti_post)
npts = 250
rr = np.linspace(rr.min(), rr.max(), npts)
write=False
if write:
    pred = np.zeros((npts, x, y, z))
    for i in range(x):
        print i, 'of ', x
        for j in range(y):
            for k in range(z):
                pred[:, i, j, k] = model.forward3(rr, Lpost[i], dpost[i], Fpost[j], Ti_post[k], mu, w, V=V_post[k])*A_post[k]
    with open("ti_solver_pred_data.p", 'wb') as outfile:
    #with open("ti_solver_pred_data_1.p", 'wb') as outfile:
        pickle.dump(pred, outfile)
else:
    #with open("ti_solver_pred_data_1.p", 'rb') as infile:
    with open("ti_solver_pred_data.p", 'rb') as infile:
        pred = pickle.load(infile)

pred = pred.reshape((npts, x*y*z))

pred_min = np.min(pred, axis=1)
pred_max = np.max(pred, axis=1)
fig, ax = plt.subplots()
ax.errorbar(r, s, fmt='.', yerr=0.03*s, color='b')
ax.fill_between(rr, pred_min, pred_max, color='r', alpha=0.5)
plt.show()
#print 'Ti (eV): ',Ti
#print 'V (km/s): ',V
#F = np.random.choice(Fpost)
#
#i = range(len(Lpost))
#j = np.random.choice(i)
#
#L = Lpost[j]
#d = dpost[j]
#
#pred = model.forward3(rr, L, d, F, Ti, mu, w, V=V)*A
#
#plt.plot(r, s, 'r')
#plt.plot(rr, pred, 'b')
#plt.show()
