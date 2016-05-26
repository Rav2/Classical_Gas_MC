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
    :return: list of lists of particle coordinates:  r=[[x1, y1, z1],...,[xN, yN, zN]]
    """
    seed()
    r = []
    for particle in range(0,N):
        xyz = [gauss(x_av, sigma), gauss(y_av, sigma), gauss(z_av, sigma)]
        r.append(xyz)
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
    energy = []
    iter = []
    for no in range(0, steps):
        k = randint(0, N-1)  #particle number
        part = r[k]
        r_cp = r[:]
        new_part = []
        for ii in range(0, 3):
            #TODO: dodac ograniczenia pudelka
            new_part.append(gauss(part[ii], sigma))
        r_cp[k]=new_part

        P_old = exp_fact(k_b, T, r, m, w)
        energy.append(exp_fact(k_b, T, r, m, w))
        iter.append(no)
        P_new = exp_fact(k_b, T, r_cp, m, w)
        P = P_new/P_old
        print(P)
        #print(P)
        if random() < P:
            r[k] = new_part
            accepted += 1
        if no > 10 and float(accepted)/float(no) < 0.5:
            sigma -= 1.0/float(no)
        elif no > 10 and float(accepted)/float(no) > 0.5:
            sigma += 1.0/float(no)
        if sigma < 0:
            sigma=0.1
    plt.plot(iter, energy)
    plt.show()
    print('acceptance: ',float(accepted) / float(steps))
    E = V(r, w, m)
    #print('roznica', E_old/E)
    print('error ', (1 - E/(1.5*N*k_b*T)) * 100, "%")






def main():
    N = 100
    m = 1
    k = 1
    T = 100
    w = 2
    sigma0 = 0.1
    steps = 10000
    r = generate_system(N, 10, 10, 10, 1)
    dynamics(r, steps, sigma0, k, T, m, w)


main()
