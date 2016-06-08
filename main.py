#!/usr/bin/python

from random import *
from fun import *
from matplotlib import pyplot as plt

def generate_system(N, x_av, y_av, z_av, sigma):
    """

    :param N: number of particles
    :param x_av: for gauss distribution
    :param y_av: for gauss distribution
    :param z_av: for gauss distribution
    :param sigma: for gauss distribution
    :return: numpy array of numpy arrays of particle coordinates:  r=[[x1, y1, z1],...,[xN, yN, zN]]
    """
    seed()
    r = np.zeros((N,3))
    for particle in range(0,N):
        xyz = np.array([gauss(x_av, sigma), gauss(y_av, sigma), gauss(z_av, sigma)])
        r[particle] = xyz
    return r


def dynamics(r, steps, sweep, sigma0, k_b, T, m, w):
    """

    :param r: vector of all particles positions
    :param steps:
    :param sigma0:
    :param k_b:
    :param T:
    :param m:
    :param w:
    :return: potential energy expected value
    """
    N = len(r)
    seed()
    sigma = sigma0
    accepted = 0
    e_factors = np.zeros((1, steps), dtype=np.float64)  # e_factors are all computed values of energy
    step_no = [x for x in range(0, steps)]  # regular list not numpy!
    cutoff = 0
    print("calculating <E>")
    for ii in range(0, steps):
        accdata = []
        for iii in range(sweep): # add energy only after sweep steps
            no = randint(0, N-1)  # particle number
            part = r[no]
            r_cp = np.array(r)
            new_part = np.zeros((1, 3))
            for jj in range(0, 3):
                new_part[0][jj] = (gauss(part[jj], sigma))
            r_cp[no] = new_part
            P = exp_fact2(k_b, T, r, r_cp, m, w)
            if random() < P:
                r = r_cp
                accepted += 1
                accdata.append(1)
            else:
                accdata.append(0)
                
        e_factors[0][ii] = V(r, w, m)
        
        # GW - better method for updating sigma
        accrate=1.0*sum(accdata)/len(accdata)
        sigma = sigma + 0.01*(accrate-0.5)
        if sigma<=0.0: sigma=1.0e-9 # cannot be negative
        #print ii, accrate, sigma

        '''
        #changing cutoff
        cutoff=0.0
        if ii >= 100 and ii % 10 == 0 and condition:
            p = np.polyfit(step_no[ii - 100:ii], e_factors[0][ii - 100:ii], 1)
            print(step_no[ii], ' ', p[0])
            if (np.sum(e_factors[0][ii - 20:ii]) - np.sum(e_factors[0][ii - 40:ii - 20])) < 100:
                cutoff = (ii) / steps  # how much beginning values we have to cut off
                print('cutoff: ', cutoff)
                condition = False'''
        cutoff = 0.1

        if ((ii+1)/steps*100) % 10 == 0:
            print(((ii+1)/steps*100), '%')

    energy = e_factors[0][int(cutoff*steps):]  # we take only energies in equilibrium state
    print("calculating u(<E>)")
    estimator = np.mean(energy)
    e_error = energy_error(energy, estimator, int((1 - cutoff) * steps))
    print('Energy error: ', e_error)
    print('last sigma: ', sigma)
    print('acceptance: ',float(accepted) / float(steps*sweep))  # should be close to 0.5
    print('Energy with expanded uncertainty(K=2): ', estimator,'+/-', 2*e_error)
    print('Theoretical value: ', 3./2. * N * k_b * T)

    plt.plot(step_no, e_factors[0])  # plotting E(steps)
    plt.ylabel("potential energy")
    plt.xlabel("MC steps")
    plt.show()


def main():
    N = 100
    m = 1.
    k = 1.
    T = 10.
    w = 0.5
    sigma0 = 1.5
    steps = 10000  # min 20
    sweep = 300
    r = generate_system(N, 0., 0., 0., 1.)
    dynamics(r, steps, sweep, sigma0, k, T, m, w)


main()
