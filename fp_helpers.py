import os
import numpy as np
from os.path import join
def make_directory(save_directory):
    if not os.path.isdir(save_directory):
        try:
            os.mkdir(save_directory)
        except OSError, e:
            print "Error making directory, maybe something else is currently making it.  Moving on gracefully..."


def peak_calculator(L, d, w, order):
    m = 2.e6 * d / w
    m0 = np.floor(m)
    #return L*np.sqrt(2.0*(1 - (m0-order)/m))
    return L * np.sqrt( m**2 / (m0 - order)**2 - 1.0)

def read_Ld_results(Ld_directory):
    fname = join(Ld_directory, "fp_Ld_post_equal_weights.dat")
    post = np.loadtxt(fname, ndmin=2)
    L = post[:, 0]
    d = post[:, 1]
    return L, d

def read_finesse_results(folder):
    fname = join(folder, "fp_F_post_equal_weights.dat")
    post = np.loadtxt(fname, ndmin=2)
    F = post[:, 0]
    Amp = post[:, 1]
    rel = post[:, 2]
    Ti_Ar = post[:, 3]
    return F, Amp, rel, Ti_Ar

def read_Ar_Ti_results(folder):
    fname = join(folder, "fp_post_equal_weights.dat")
    post = np.loadtxt(fname, ndmin=2)
    Ti = post[:, 0]
    A = post[:, 1]
    V = post[:, 2]
    return Ti, A, V





