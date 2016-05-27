from math import floor, exp, ceil
import numpy as np
def V(r, w, m):
    """
    :param r: list of position vectors of all particles
    :param w: omega of oscillator
    :param m:  mass of particles(all have equal one)
    :return: potential energy
    """
    return 0.5 * m * (w ** 2)*np.sum(np.square(r))


def exp_fact(k, T, r, m, w):
    """
    :param k: Boltzmann's const
    :param T: temperature
    :param r: list of position vectors of all particles
    :param w: omega of oscillator
    :param m:  mass of particles(all have equal one)
    :return: exponential factor
    """
    beta = 1.0/(k*T)
    U = V(r, w, m)
    val = exp(-beta*U)
    return val

