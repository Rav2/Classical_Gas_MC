from math import exp
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


def energy_error(energies, left_steps):
    """
    :param energies: array of energy factors taken at equillibrium state
    :param left_steps: number of particles
    :return: delta_E energy error
    """

    estimator_E = np.mean(energies)
    R_0 = (np.std(energies))**2
    k_range = np.arange(1, left_steps)
    R_k = np.zeros_like(k_range, dtype=float)
    for k in k_range:
        if k % int(0.1 * left_steps)==0:
            print(k / left_steps * 100 + 10, '%')
        #usunac +1
        R_k[k-1] = np.sum((energies[0:left_steps - k] - estimator_E) * (energies[k:left_steps+1] - estimator_E))
        #R_k[k-1] = np.sum((energies[k_range < left_steps - k] - estimator_E) * (energies[k_range > k] - estimator_E))
    R_k = R_k[R_k > 0.1]

    delta_E = np.sqrt(R_0 * (1 + 2 * sum(R_k)) / left_steps)
    return delta_E

