import numpy as np
import scipy.special


def UniformPrior(cube, a, b):
    """Maps Uniform distribution [0, 1] to a Uniform distribution [a, b]

    :param cube: number in the range [0, 1]
    :param a: left edge of Uniform distribution
    :param b: right edge of Uniform distribution
    :returns: Maps cube in [0,1] to [a,b]
    """
    return a + (b-a) * cube


def LogUniformPrior(cube, a, b):
    """Maps Uniform distribution [0, 1] to a Uniform distribution in log10 space from log10(a) to log10(b)

    :param cube: number in the range [0, 1]
    :param a: left edge of LogUniform distribution in real space
    :param b: right edge of LogUniform distribution in real space
    :returns: Maps cube in [0,1] to [a,b] with a LogUniform distribution
    """

    return 10.0**UniformPrior(cube, a, b)


def NormalPrior(cube, mu, sigma):
    """Maps [0,1] to normal distribution with mean mu and variace sigma**2

    :param cube: number in the range [0, 1]
    :param mu: mean of the normal distribution
    :param sigma: standard deviation of the normal distribution (variance=sigma**2)
    :returns: Maps [0,1] to the normal distribution
    """
    val = np.sqrt(2.0) * scipy.special.erfinv(2*cube-1)
    return sigma*val+mu 


def NormalCDF(x):
    """Cumulative distribution function for the Normal distribution

    :param x: number between 0 and 1
    :returns: the mapping value from the unit normal distribution
    """
    return 0.5 * (1.0 + scipy.special.erf(x/np.sqrt(2.0)))


def TruncatedNormal(cube, mu, sigma, a, b):
    """Maps [0,1] to the truncated normal distribution with mean mu, variace sigma**2, left edge a, and right edge b

    :param cube: number in the range [0, 1]
    :param mu: mean of the normal distribution
    :param sigma: standard deviation of the normal distribution (variance=sigma**2)
    :param a: left edge for truncation
    :param b: right edge for truncation
    :returns: Maps [0,1] to the truncated normal distribution
    """

    alpha = (a-mu)/sigma
    beta = (b-mu)/sigma

    Phi_alpha = NormalCDF(alpha)
    Phi_beta = NormalCDF(beta)
    Z = Phi_beta - Phi_alpha

    val = scipy.special.erfinv(2.0*Z*cube + 2*Phi_alpha -1.0)

    return np.sqrt(2.0) * sigma * val + mu


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    cubes = np.random.rand(10000)
    vals = TruncatedNormal(cubes, 0.0, 0.75, -1000, 0.5)

    fig, ax = plt.subplots()
    ax.hist(vals, bins=30)
    plt.show()



