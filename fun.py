from math import exp


def T(p, m):
    """
    :param p: list of momentum vectors of all particles
    :param m: mass of particles(all have equal one)
    :return: kinetic energy of system
    """
    e = 0.0
    for k in p:
        for p_i in k:
            e += p_i**2/(2*m)
    return e


def V(r, w, m):
    """
    :param r: list of position vectors of all particles
    :param w: omega of oscillator
    :param m:  mass of particles(all have equal one)
    :return: potential energy
    """
    e = 0.0
    for x in r:
        for x_i in x:
            e += x_i ** 2 / (2 * m)
    return e


def H( p, r, m, w):
    """
    :param p: list of momentum vectors of all particles
    :param r: list of position vectors of all particles
    :param w: omega of oscillator
    :param m:  mass of particles(all have equal one)
    :return: total energy (Hamiltonian)
    """
    return T(p,m)+V(r,w,m)


def exp_fact(k, T, p, r, m, w):
    """
    :param k: Boltzmann's const
    :param T: temperature
    :param p: list of momentum vectors of all particles
    :param r: list of position vectors of all particles
    :param w: omega of oscillator
    :param m:  mass of particles(all have equal one)
    :return: exponential factor
    """
    beta = 1/(k*T)
    return exp(-beta*H(p, r, m, w))

