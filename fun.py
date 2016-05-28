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


def energy_error(energies, left_steps):
    """
    :param energies: array of energy factors taken at equillibrium state
    :param left_steps: number of particles
    :return: delta_E energy error
    """

    estimator_E = np.mean(energies)
    R_0 = (np.std(energies))**2
    R_k = np.array([])
    for k in range(1, left_steps):
        if k % int(0.1 * left_steps)==0:
            print(k / left_steps * 100 + 10, '%')
        rk = 0

        for i in range(1, left_steps-k):
            rk += (energies[i] - estimator_E) * (energies[i+k] - estimator_E)
        #rk = sum()
        if rk > 0.1:
            R_k = np.append(R_k,rk)

    delta_E = np.sqrt(R_0 * (1 + 2 * sum(R_k)) / left_steps)
    return delta_E