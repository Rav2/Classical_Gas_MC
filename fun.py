from math import exp
import numpy as np
def V(r, w, m):
    """
    :param r: list of position vectors of all particles
    :param w: omega of oscillator
    :param m:  mass of particles(all have equal one)
    :return: potential energy
    """
    return 0.5 * m * w * w * np.sum(np.square(r))

def V_len(r, eps, sigma, w, m):
    """
    Calculates potential energy of system of particles for  Lennard-Jones potential field
    :param r: list of position vectors of all particles
    :param eps: parameter of Lennard-Jones potential
    :param sigma: 2nd parameter of Lennard-Jones potential
    :return: potential energy for all particles
    """
    sum = 0.0
    for part in r:
        sum += sigma**12/((np.sum(part*part))**6) - sigma**6/((np.sum(part*part))**3)
    return 4*eps*sum +  V(r, w, m)

def exp_fact2_len(eps, sigma, rold, rnew, k, T, w, m):
    """
    Calculates value of the exponential function for two states of system
    :param eps: parameter of Lennard-Jones potential
    :param sigma: 2nd parameter of Lennard-Jones potential
    :param rold: particle position before change
    :param rnew: particle position after change
    :param k: Boltzmann's const
    :param T: temperature
    :return: exponential factor for Lennard-Jones potential field
    """

    beta = 1.0/(k*T)
    U = V_len(rnew, eps, sigma, w, m) - V_len(rold, eps, sigma,w, m)
    # print(-beta*U)
    try:
        val = exp(-beta*U)
    except OverflowError:
        val = float('inf')
    return val

def exp_fact2(k, T, rold, rnew, m, w):
    """
    :param k: Boltzmann's const
    :param T: temperature
    :param r: list of position vectors of all particles
    :param w: omega of oscillator
    :param m:  mass of particles(all have equal one)
    :return: exponential factor
    """
    beta = 1.0/(k*T)
    U = V(rnew, w, m)-V(rold, w, m)
    val = exp(-beta*U)
    return val


def energy_error(energies, estimator_E, left_steps):
    """
    :param energies: array of energy factors taken at equillibrium state
    :param left_steps: number of particles
    :return: delta_E energy error
    """
    R_0 = np.var(energies)
    k_range = np.arange(1, left_steps)
    R_k = np.zeros_like(k_range, dtype=float)
    energies = energies - estimator_E
    for k in k_range:
        if ((k+1)/left_steps*100)%20 == 0:
            print((k+1)/left_steps*100, '%')
        R_k[k-1] = np.sum((energies[0:left_steps - k]) * (energies[k:left_steps+1]))/((left_steps)) / R_0
    R_k = R_k[R_k > 0.1]

    delta_E = np.sqrt(R_0*(1.0 + 2 * sum(R_k)) / left_steps)
    return delta_E

