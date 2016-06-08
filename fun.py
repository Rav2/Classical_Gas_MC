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
        if (k/left_steps*1000)%10 == 0:
            print(k/left_steps*100, '%')
        # if k % int(0.1 * left_steps)==0:
        #     print(k / left_steps * 100 + 10, '%')
        #usunac +1
        R_k[k-1] = np.sum((energies[0:left_steps - k]) * (energies[k:left_steps+1]))/((left_steps+1)-k)
        #R_k[k-1] = np.sum((energies[k_range < left_steps - k] - estimator_E) * (energies[k_range > k] - estimator_E))
    R_k = R_k[R_k > 0.1]

    delta_E = np.sqrt((R_0 + 2 * sum(R_k)) / left_steps)
    return delta_E

