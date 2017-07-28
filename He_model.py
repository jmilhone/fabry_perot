from __future__ import division
from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'
import model2

wavelengths = [468.5376849,
               468.54072254,
               468.55244041,
               468.55680062,
               468.57038494,
               468.5704380,
               468.575707958,
               468.5757974,
               468.58040922,
               468.58308897 ,
               468.588412282,
               468.59055531,
               468.591788438]

#wavelengths = [1945.3, 1807.1, 1273.5, 1075.0, 456.0, 453.9, 214.1, 209.0, 0.0, -122.0, -364.4, -461.9, 518.1]
#wavelengths = [x*.001 + 21335.0822 for x in wavelengths]
#wavelengths = [x*100.0 for x in wavelengths]
#wavelengths = [1.0 / x for x in wavelengths]
#wavelengths = [x * 1e9 for x in wavelengths]

#for x in wavelengths:
#    print x

Aij = [9.3894e+07,
       4.9071e+07,
       9.7947e+06,
       4.9070e+07,
       2.0600e+08,
       1.1266e+08,
       5.5631e+05,
       1.8776e+07,
       2.2070e+08,
       1.4713e+07,
       5.0067e+06,
       1.9588e+07,
       5.5636e+06]


#groupings = [[2,11],
#             [1,3,6,10,12],
#             [0,5,7],
#             [4,8,9]]
#
##            0  1  2  3  4  5  6  7  8  9  10 11 12
#map2group = [2, 1, 0, 1, 3, 2, 1, 2, 3, 3, 1, 0, 1]  

groupings = [[0, 5, 7],
             [1, 3],
             [2,11],
             [4, 8, 9],
             [6, 10, 12]]
map2group = [0, 1, 2, 1, 3, 0, 4, 0, 3, 3, 4, 2, 4]

def calculate_branching_fractions():
    branching_norms = []
    for group in groupings:
        branching_norms.append(sum(Aij[x] for x in  group))

    branching_fractions = range(13)
    for idx, group in enumerate(groupings):
        for w in group:
            branching_fractions[w] = Aij[w] / branching_norms[idx]

    return branching_fractions


def main(L, d, V=0, Ti=0.5):
    bf = calculate_branching_fractions()
    amplitudes = [0.32, 0.2, 1.13, 1.61, 1.0]
    #amplitudes = [0.32, 0.2, 1.13, 1.61]
    individual_amplitudes = [bf[x] * amplitudes[map2group[x]] for x in range(13)]
    print bf
    print individual_amplitudes

    #w = np.linspace(468.5, 468.625, 10000)
    w = np.linspace(min(wavelengths)-.01, max(wavelengths)+.01, 10000)
    filt = model2.eval_spec(w, 1.0, 468.67, .46/2.0 /np.sqrt(2.0 * np.log(2.0)))
    filt *= .3423 / filt.max()
    spec = None
    #Ti = 0.46
    mu = 4.0
    rarr = np.linspace(0.0, 2000.0, 4000)
    

    lin_out = model2.forward5(rarr, L, d, 22.7, Ti*np.ones(13), 
            4.0*np.ones(13), individual_amplitudes, wavelengths, V=V)
    #plt.plot(rarr**2, lin_out)
    #plt.show(block=False)
    #for w0, amp in zip(wavelengths, individual_amplitudes):
    #    sigma = 3.276569e-5 * np.sqrt(Ti/mu) * w0
    #    temp = model2.eval_spec(w, amp, w0, sigma, V=0)
    #    #temp *= filt
    #    if spec is None:
    #        spec = temp
    #    else:
    #        spec = np.vstack((spec, temp))
    #tot_spec = np.sum(spec, axis=0)
    #fig, ax = plt.subplots()
    #ax.plot(w, tot_spec)
    #for i in range(12):
    #    ax.plot(w, spec[i, :])
    #plt.show()

    return rarr, lin_out

if __name__ == "__main__":
    L = 149.496800727 / .004
    d = 0.88342569509 
    rarr, lin_out = main(L, d)
    #plt.plot(rarr**2, lin_out)
    #plt.show(block=False)
    #for w0, amp in zip(wavelengths, individual_amplitudes):
    #    sigma = 3.276569e-5 * np.sqrt(Ti/mu) * w0
    #    temp = model.eval_spec(w, amp, w0, sigma, V=0)
    #    temp *= filt
    #    if spec is None:
    #        spec = temp
    #    else:
    #        spec = np.vstack((spec, temp))
    #tot_spec = np.sum(spec, axis=0)
    #fig, ax = plt.subplots()
    #ax.plot(w, tot_spec)
    #for i in range(12):
    #    ax.plot(w, spec[i, :])
    #plt.show()


