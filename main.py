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


def dynamics(r, steps, sigma0, k_b, T, m, w):
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
    condition = True
    cutoff = 0
    for ii in range(0, steps):
        no = randint(0, N-1)  # particle number
        part = r[no]
        r_cp = np.array(r)
        new_part = np.zeros((1, 3))
        for jj in range(0, 3):
            new_part[0][jj]=(gauss(part[jj], sigma))
        r_cp[no] = new_part
        P_old = exp_fact(k_b, T, r, m, w)
        P_new = exp_fact(k_b, T, r_cp, m, w)
        P = P_new/P_old
        # print(P)
        if random() < P:
            r = r_cp
            accepted += 1
        e_factors[0][ii] = V(r, m, w)

        #changing cutoff
        if ii >= 40 and ii % 10 == 0 and condition:
            if (np.sum(e_factors[0][ii - 20:ii]) - np.sum(e_factors[0][ii - 40:ii - 20])) < 100:
                #print(np.sum(e_factors[0][ii - 20:ii]) - np.sum(e_factors[0][ii - 40:ii - 20]))
                cutoff = (ii) / steps  # how much beginning values we have to cut off
                print('cutoff: ', cutoff)
                condition = False

        #changing the value of sigma dynamically
        if(ii < 1000):
            kk=ii
        else:
            kk=1000
        if ii > 10 and accepted/ii < 0.5:
            sigma -= 1.0/kk
        elif ii > 10 and accepted/ii > 0.5:
            sigma += 1.0/kk
        if sigma < 0:
            sigma = 0.1
        #printing progress
        if ((ii+1)/steps*100) % 100 == 0:
            print(((ii+1)/steps*10), '%')
    energy = e_factors[0][int(cutoff*steps):]  # we take only energies in equilibrium state

    e_error = energy_error(energy, int((1 - cutoff) * steps))
    print('Energy error: ', e_error)
    print('last sigma: ', sigma)
    print('acceptance: ',float(accepted) / float(steps))  # should be close to 0.5
    print('Energy with expanded uncertainty(K=2): ',np.mean(energy),'+/-', 2*e_error)

    plt.plot(step_no, e_factors[0])  # plotting E(steps)
    plt.ylabel("potential energy")
    plt.xlabel("MC steps")
    plt.show()

    #E = sum(energy)/(int(cutoff*steps))  # mean value of energy as estimator
    #print('error: ', (1 - E/(1.5*N*k_b*T)) * 100, "%")  # relative error, only for us



def main():
    N = 10000
    m = 1
    k = 1
    T = 100 # increase when something is not working as it should
    w = 1.5
    sigma0 = 1
    steps = 10000  # min 20
    r = generate_system(N, 0, 0, 0, 1)
    dynamics(r, steps, sigma0, k, T, m, w)


main()
