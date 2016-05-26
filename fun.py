from math import exp



def V(r, w, m):
    """
    :param r: list of position vectors of all particles
    :param w: omega of oscillator
    :param m:  mass of particles(all have equal one)
    :return: potential energy
    """
    E = 0.0
    for q in r:
        for x_i in range(0, len(q)):
            E += q[x_i]**2
    E *= 0.5 * m * (w**2)
    return E


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
    val = exp(-beta*V(r, m, w))
    return val

