import matplotlib.pyplot as plt
import numpy as np
from forward_model import main

_, arr_dict = main('Ar', plotit=False)
print 'done with one'
_, arr_dict_sa = main('Ar', plotit=False, sm_ang=True)
print 'plotting'
f, ax = plt.subplots(2, 1, figsize=(16, 9))
ax[0].plot(arr_dict['pixel_arr'], arr_dict['lin_data'], color='r', label='full')
ax[0].plot(arr_dict_sa['pixel_arr'], arr_dict_sa['lin_data'], color='b', label='sm_angle')
ax[0].plot([2010]*2, [0, 1.05], 'g--', label='edge of ccd')
ax[0].set_ylim([0, 1.05])
ax[0].legend(loc='best')
ax[1].plot(arr_dict['pixel_arr'], arr_dict['lin_data']-arr_dict_sa['lin_data'])
ax[1].plot([2010]*2, [-1, 1], 'g--', label='edge of ccd')
ax[1].set_ylim([-1, 1])
plt.show(block=False)

x = np.linspace(0.,7.8,1000)
full = 1./np.sqrt(1+(x/350.)**2)
appr = 1.-0.5*(x/350.)**2
fact34 = 2.*0.88e6* (1./3400)
fact40 = 2.*0.88e6* (1./4000)
f,ax = plt.subplots(figsize=(16,9))
ax.fill_between(x, fact40*(full-appr)*1.e6, fact34*(full-appr)*1.e6, color='b', alpha=0.7)
ax.set_xlabel('r (mm)')
ax.set_ylabel(r'$\Delta\lambda$ (fm)')
plt.show()
