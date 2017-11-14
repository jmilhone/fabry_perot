import matplotlib.pyplot as plt
import numpy as np
from forward_model import main

# _, arr_dict = main('Ar', plotit=False)
# print 'done with one'
# _, arr_dict_sa = main('Ar', plotit=False, sm_ang=True)
# print 'plotting'
# f, ax = plt.subplots(2, 1, figsize=(16, 9))
# ax[0].plot(arr_dict['pixel_arr'], arr_dict['lin_data'], color='r', label='full')
# ax[0].plot(arr_dict_sa['pixel_arr'], arr_dict_sa['lin_data'], color='b', label='sm_angle')
# ax[0].plot([2010]*2, [0, 1.05], 'g--', label='edge of ccd')
# ax[0].set_ylim([0, 1.05])
# ax[0].legend(loc='best')
# ax[1].plot(arr_dict['pixel_arr'], arr_dict['lin_data']-arr_dict_sa['lin_data'])
# ax[1].plot([2010]*2, [-1, 1], 'g--', label='edge of ccd')
# ax[1].set_ylim([-1, 1])
# plt.show(block=False)

d_arr = np.array([0.87,0.88,0.89])
f_arr = np.array([148.,150.,152.])
lam_arr = np.array([468.6,488.])
j_arr = np.array([0,1,2,3,4])
d_lam = 1.e-6 * np.linspace(0., 0.1, 1000)
f, ax = plt.subplots(figsize=(16, 9))
for d in d_arr:
    for f in f_arr:
        for lam in lam_arr:
            for j in j_arr:
                m0 = np.floor(2.*d*1e6 / lam)
                m0j = m0-j
                r_full = f * np.sqrt(((2.*d*1e6 / lam) / m0j)**2 - 1.)
                r_appr = f * np.sqrt(2 - ((m0j*lam)/(1e6*d)))
                full = r_full - f*np.sqrt(((2*d*np.sqrt(f**2+r_full**2))/(2*d*f+d_lam*m0j*np.sqrt(f**2+r_full**2)))**2 - 1.)
                appr = r_appr - np.sqrt(r_appr**2 -(1./d)*m0j*d_lam*f**2)
                ax.fill_between(d_lam*1.e6, full, appr, lw=2, label='d={0}, f={1}, lam={2}, j={3}'.format(d, f, lam, j))
ax.set_xlabel(r'$\Delta\lambda$ (nm)')
ax.set_ylabel('Full - Approx (mm)')
#ax.legend(loc='best', fontsize='8')
plt.show()