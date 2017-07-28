from __future__ import division
import numpy as np
#import model2
import model
# vacuum
#wavelengths = [1945.3, 1807.1, 1273.5, 1075.0, 456.0, 453.9, 214.1, 209.0, 0.0, -122.0, -364.4, -461.9, 518.1]
#wavelengths = [x*.001 + 21335.0822 for x in wavelengths]
#wavelengths = [x*100.0 for x in wavelengths]
#wavelengths = [1.0 / x for x in wavelengths]
#wavelengths = [x * 1e9 for x in wavelengths]

#for x in wavelengths:
#    print x

# air
wavelengths = np.array([468.5376849,
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
               468.591788438])

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

def peak_locations(L, d, w, orders=1):
    m = 2.0e6*d / w
    m0 = np.floor(m)

    if orders == 1:
        return L*np.sqrt(m**2 / m0**2 - 1.0)
    else:
        return [L*np.sqrt(m**2 / (m0-j)**2 - 1.0) for j in range(orders)]

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

bf = calculate_branching_fractions()


def model_output(r, L, d, Ti, amplitudes, V=0.0, mu=4.0026, F=22.7):
    n = 13
    individual_amplitudes = [bf[x]*amplitudes[map2group[x]] for x in range(n)]
    # going to attempt to normalize the amplitudes 
    #tot = np.sum(individual_amplitudes)
    #individual_amplitudes = [x / tot for x in individual_amplitudes] 
    #lin_out = model2.forward5(r, L, d, F, Ti*np.ones(n), mu*np.ones(n), 
    #        individual_amplitudes, wavelengths, V=V)
    lin_out = model.forward4(r, L, d, F, Ti*np.ones(n), mu*np.ones(n), 
            individual_amplitudes, wavelengths.copy(), V=V)
    return lin_out






