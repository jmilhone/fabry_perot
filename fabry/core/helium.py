from __future__ import division, print_function
import numpy as np

wavelength = [468.5376849, 468.54072254, 468.55244041, 468.55680062,
              468.57038494, 468.5704380, 468.575707958, 468.5757974,
              468.58040922, 468.58308897, 468.588412282, 468.59055531,
              468.591788438]

upper_states = ['4d', '4p', '4s', '4p',
                '4f', '4d', '4p', '4d',
                '4f', '4f', '4p', '4s',
                '4p']

Aij = [9.3894E7, 4.9071e7, 9.7947e6, 4.9070e7, 2.0600e8, 1.1266e8, 5.5631e5,
       1.8776e7, 2.2070e8, 1.4713e7, 5.0067e6, 1.9588e7, 5.5636e6]

upper_to_lower_Aij = {'4s': np.array([1.3758e7, 2.7515e7, 9.7947e6, 1.9588e7]),
                      '4p': np.array([1.5478e8, 1.5478e8, 1.0915e9, 1.0915e9, 4.9071e7, 4.907e7, 5.5631e5, 5.0067e6,
                                      5.5636e6]),
                      '4d': np.array([2.7518e8, 3.3015e8, 5.5024e7, 9.3894e7, 1.1266e8, 1.8776e7]),
                      '4f': np.array([2.06e8, 2.207e8, 1.4713e7, ]),
                      }


def branching_fractions():
    gamma = np.zeros_like(wavelength)

    totals = {key: np.sum(einstein_coefficients) for (key, einstein_coefficients) in upper_to_lower_Aij.items()}

    for idx, (upper, einstein_coefficient) in enumerate(zip(upper_states, Aij)):
        gamma[idx] = einstein_coefficient / totals[upper]

    gamma = np.array(gamma)
    gamma /= gamma.max()

    return gamma
