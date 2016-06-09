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

def V_len(r, eps, sigma):
    """

    :param r:
    :param eps:
    :param sigma:
    :return:
    """
    sum = 0.0
    for part in r:
        sum += sigma**12/((part[0]**2+part[1]**2+part[2]*2)**6) - sigma**6//((part[0]**2+part[1]**2+part[2]*2)**3)
    return 4*eps*sum

def exp_fact2_len(eps, sigma, rold, rnew, k, T):

    beta = 1.0/(k*T)
    U = V_len(rnew, eps, sigma)-V_len(rold, eps, sigma)
    # print(U)
    val = exp(-beta*U)
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
        if (k/left_steps*100)%10 == 0:
            print(k/left_steps*100, '%')
        R_k[k-1] = np.sum((energies[0:left_steps - k]) * (energies[k:left_steps+1]))/((left_steps)) / R_0
    R_k = R_k[R_k > 0.1]

    delta_E = np.sqrt(R_0*(1.0 + 2 * sum(R_k)) / left_steps)
    return delta_E

