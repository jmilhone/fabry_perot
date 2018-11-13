from __future__ import division, print_function
#import sys
#sys.path.append("../")
import numpy as np
import matplotlib.pyplot as plt
from fabry.plasma import plasma


def main():
    impact_factor = 0.3
    Ti = 0.5
    w0 = 488.0
    mu = 40.0
    Lnu = 50.0
    Vouter = 10.0
    rmax = 40.0
    nr = 201
    nlambda=2000
    Lne = 2.0
    R_outer = 35.0

    r = np.linspace(0.0, rmax+5, nr)
    v = plasma.pcx_velocity_profile(r, Lnu, R_outer, Vouter)

    w, spec = plasma.calculate_pcx_chord_emission(impact_factor, Ti, w0, mu, Lnu, Vouter, nr=401)

    fig, ax = plt.subplots()
    ax.plot(w, spec)
    plt.show()

if __name__ == "__main__":
    main()
