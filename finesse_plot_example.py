from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import model
import plottingtools.core as ptools

rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'



L = 150.0 / .004
d = 0.88

w = 488.0

Farr = [2, 5, 15]

r = np.linspace(0.0, 1200.0, 10000)
cos_theta = L / np.sqrt(L**2 + r**2)

fig, ax = plt.subplots()
for F in Farr:
    temp = model.eval_airy(w, cos_theta, d, F)
    ax.plot(cos_theta, temp, label=str(F))

lg = ax.legend(frameon=False, fontsize=18, ncol=3, title="Finesse")
lg.get_title().set_fontsize(18)
ax.set_xlabel(r"$\cos{\theta}$", fontsize=20)
ax.set_ylabel("Amplitude (A.U.)", fontsize=20)
ptools.add_thick_box(ax, minor=False)
ax.set_ylim(0.0, 1.25)
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.tight_layout()
fname = "Plots/FinesseExample.pdf"
fig.savefig(fname)
plt.show()

