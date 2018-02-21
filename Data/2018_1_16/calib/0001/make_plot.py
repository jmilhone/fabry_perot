import numpy as np
import matplotlib.pyplot as plt
from tools.file_io import read_Ld_results, h5_2_dict
from core.models import peak_calculator
from tools.plotting import my_hist,peak_plot

bins = 30

L_real = 149.883
d_real = 0.881235
px = 0.004
wavelengths = {'Th':487.873302, 'Ar':487.98634}

data = h5_2_dict('input_data.h5')
Lpost, dpost = read_Ld_results('.')

hists = {}
r_pks = {}
for w, pk in data['peaks'].items():
    h_list = []
    r_pk = []
    for i,n in enumerate(data['orders'][w]):
        h_list.append(peak_calculator(Lpost,dpost,float(w),n))
        r_pk.append(peak_calculator(L_real/px, d_real, float(w), n))
    hists[w] = h_list
    r_pks[w] = r_pk

norder = max([len(x) for x in data['orders'].values()])
nwaves = len(data['peaks'].keys())

fig,ax = plt.subplots(norder, nwaves, figsize=(12,10))
axx = ax.reshape((norder,nwaves))
for i,w in enumerate(hists.keys()):
    for j,hist in enumerate(hists[w]):
        my_hist(axx[j,i],hist,bins=bins)
        axx[j,i].axvline(data['peaks'][w][j], color='k')
        axx[j,i].axvspan(data['peaks'][w][j] - data['pkerr'][w][j], data['peaks'][w][j] + data['pkerr'][w][j], color='k', alpha=0.2)
        axx[j,i].axvline(r_pks[w][j], color='r')
        axx[j,i].set_ylabel('Order {0:d}'.format(data['orders'][w][j]))
    axx[0,i].set_title(data['names'][w],fontsize=18)
    axx[-1,i].set_xlabel('R (px)')
for ax in axx.flatten():
    ax.get_xaxis().get_major_formatter().set_scientific(False)
plt.show(block=False)


fig1,ax1 = plt.subplots(2,1,figsize=(6,8))
my_hist(ax1[0],Lpost*0.004,bins=bins)
ax1[0].set_xlabel('L (mm)')
ax1[0].set_ylabel('P (L)')
ax1[0].axvline(L_real, color='r')
my_hist(ax1[1],dpost,bins=bins)
ax1[1].set_xlabel('d (mm)')
ax1[1].set_ylabel('P (d)')
ax1[1].axvline(d_real,color='r')
for ax in ax1.flatten():
    ax.get_xaxis().get_major_formatter().set_scientific(False)

plt.show(block=False)

fig2,ax2 = plt.subplots(figsize=(10,6))
peak_plot(data['r'],data['sig'],data['peaks'],data['orders'],data['names'],fax=(fig2,ax2))
for w, pks in r_pks.items():
    for pk in pks:
        ax2.axvline(pk,color='r',linestyle='--')
plt.show()
