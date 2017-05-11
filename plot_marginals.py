from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
import plottingtools.core as ptools
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from os.path import join

params = ['L', 'd', 'Finesse', 'Ti_Th', 'Ti_Ar', 'Amp_Th', 'Amp_Ar', 'r0']
labels = ["L (mm)", "d (mm)", "Finesse", r"$T_{i, Th}$ (K)", r"$T_{i, Ar}$ (eV)", r"$A_{Th}$ (Counts)", r"$A_{Ar}$ (Counts)", r"$r_0$ (px)"]
prob_labels = ["P(L)", "P(d)", "P(Finesse)", r"P($T_{i,Th}$)", r"P($T_{i,Ar}$)", r"P($A_{Th}$)", r"P($A_{Ar}$)", r"P($r_0$)"]

#params = ["Ti", "V", "r0", "A"]
#labels = ["Ti (eV)", "V (km/s)", "r0 (px)", "A (Counts)"]
#prob_labels = ["P(Ti)", "P(V)", "P(r0)", "P(A)"]

nparams = len(params)
#fname = "saves/Ar_solver_run4/fp_full_post_equal_weights.dat"
fname = "saves/full_solver_run16/fp_full_post_equal_weights.dat"
post = np.loadtxt(fname, ndmin=2)
#post[:, 0] *= 0.004  # convert px to mm
#post[:, 3] *= 300.0 / .025 # convert to K from eV

#folder = "Plots/run15/"
save = False
for idx, param in enumerate(params):
    fig, ax = plt.subplots()
    hist, bins = np.histogram(post[:, idx], density=True, bins=40)
    bw = bins[1]-bins[0]
    ax.bar(bins[0:-1], hist*bw, width=bw)
    print np.sum(hist*bw), bw, np.sum(hist)
    #ax.plot(bins[:-1], hist)
    #ax.hist(post[:, idx], bins=bins, normed=True)
    ptools.add_labels(ax, labels[idx], prob_labels[idx])
    ptools.add_thick_box(ax, minor=False)
    print param,np.mean(post[:, idx]), np.std(post[:, idx])

    plt.tight_layout()
    if save:
        fname = join(folder, "{}_marginal.pdf".format(param))
        fig.savefig(fname)
        plt.close(fig)
    else:
        plt.show()

Lval = np.mean(post[:,0])
i = np.argmin(np.abs(Lval - post[:,0]))
print post[i, :]

for i, param1 in enumerate(params):
    for j, param2 in enumerate(params[0:i]):
        fig, ax = plt.subplots()
        hist, xx, yy = np.histogram2d(post[:, j], post[:, i], bins=100, normed=True)
        #_, _, _, im = ax.hist2d(post[:, j], post[:, i], bins=100, normed=True)
        #_, _, _, im = ax.hist2d(post[:, j], post[:, i], bins=[xx,yy], normed=True)
        dx = xx[1]-xx[0]
        dy = yy[1]-yy[0]
        im = ax.contourf(yy[0:-1], xx[0:-1], hist*dx*dy)
        ax_divider = make_axes_locatable(ax)
        cax = ax_divider.append_axes("right", size="7%", pad="2%")
        cb = plt.colorbar(im, cax=cax)
        #cax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        cax.tick_params(labelsize=18)
        #ptools.add_thick_box(cax)
        ptools.add_thick_box(ax, minor=False)
        #ptools.add_labels(ax, labels[j], labels[i])
        ptools.add_labels(ax, labels[i], labels[j])
        plt.tight_layout()
        if save:
            fname = join(folder, "{}_{}_marginal.pdf".format(param2, param1))
            fig.savefig(fname)
            plt.close(fig)
        else:
            plt.show()
#fig = plt.figure(figsize=(12,12))
#s = 1
#for i in range(nparams):
#    #ax = plt.subplot(nparams, nparams, nparams*i + i +1)
#    ax = plt.subplot(6, 6, s)
#    s += 1
#    ax.hist(post[:, i], bins=40, normed=True)
#    for j in range(i):
#        #ax = plt.subplot(nparams, nparams, nparams * j + i + 1)
#        ax = plt.subplot(6, 6, s)
#        s +=1
#        _, _, _, im = ax.hist2d(post[:, j], post[:, i], bins=100, normed=True)
#        ax_divider = make_axes_locatable(ax)
#        cax = ax_divider.append_axes("right", size="7%", pad="2%")
#        cb = plt.colorbar(im, cax=cax)
#plt.tight_layout()
#plt.show()
