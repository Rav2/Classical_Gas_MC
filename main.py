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


def dynamics(r, steps, sweep, sigma0, k_b, T, m, w, LJ_pot, LJ_eps, LJ_sigma):
    """

    :param r: vector of all particles positions
    :param steps: how many values of total energy will we save
    :param sweep: how many moves will we do in single step
    :param sigma0:
    :param k_b:
    :param T:
    :param m:
    :param w:
    :param LJ_pot: if True, Lennard-Jones potential will be used, otherwise harmonic one
    :param LJ_eps:
    :param LJ_sigma:
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
            if LJ_pot:
                P = exp_fact2_len(LJ_eps, LJ_sigma, r, r_cp, k_b, T, w, m)
            else:
                P = exp_fact2(k_b, T, r, r_cp, m, w)
            # print("P=",P)
            if random() < P:
                r = r_cp
                accepted += 1
                accdata.append(1)
            else:
                accdata.append(0)

        if LJ_pot:
            e_factors[0][ii] = V_len(r, LJ_eps, LJ_sigma, w, m)
        else:
            e_factors[0][ii] = V(r, w, m)


        # if iii > 10 and accepted / iii < 0.5:
        #     sigma -= 1.0 / iii
        # elif ii > 10 and accepted / iii > 0.5:
        #     sigma += 1.0 / iii
        # if sigma < 0:
        #     sigma = 1.e-9
        # GW - better method for updating sigma
        accrate=1.0*sum(accdata)/len(accdata)
        sigma = sigma + 0.01*(accrate-0.5)
        if sigma <= 0.0: sigma = 1.0e-9 # must be positive

        cutoff = 0.1

        if ((ii+1)/steps*100) % 10 == 0:
            print(((ii+1)/steps*100), '%')

    energy = e_factors[0][int(cutoff*steps):]  # we take only energies in equilibrium state
    print("calculating u(<E>)")
    estimator = np.mean(energy)
    e_error = energy_error(energy, estimator, len(energy))
    print("parameters: w=", w, " m=", m, " sigma0=", sigma0, " k_B=", k_b, " T=", T, " LJ=", LJ_pot, "LJ_eps=", LJ_eps, "LJ_sigma=", LJ_sigma)
    print('Energy error: ', e_error)
    print('last sigma: ', sigma)
    print('acceptance: ',float(accepted) / float(steps*sweep))  # should be close to 0.5
    print('Energy with expanded uncertainty(K=2): ', estimator,'+/-', 2*e_error)
    if not LJ_pot:
        print('Theoretical value: ', 3./2. * N * k_b * T)

    # plotting U(steps)
    # plt.plot(step_no, e_factors[0])  # plotting E(steps)
    # plt.ylabel("potential energy")
    # plt.xlabel("MC steps")
    # plt.savefig("N_",N,"T_",T,".png")
    return estimator, e_error



def main():
    N = [40, 200, 600]
    leg_entries = ['N = 40', 'N = 200', 'N = 600']
    m = 1.
    k = 1.
    T = [1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]
    w = 2.5
    sigma0 = 1.5
    lj_eps = 1.
    lj_sigma = 1.
    steps = 1000  # min 20
    sweep = 300
    #calculations and plot for harmonic potential
    # results_E = []
    # results_uE = []
    # for n in N:
    #     E = []
    #     uE = []
    #     for t in T:
    #         r = generate_system(n, 0., 0., 0., 1.)
    #         e, ue = dynamics(r, steps, sweep, sigma0, k, t, m, w, False, lj_eps, lj_sigma)
    #         E.append(e)
    #         uE.append(ue)
    #     results_E.append(E)
    #     results_uE.append(uE)
    # for ii in range(0, len(results_E)):
    #     plt.errorbar(T, 3./2.*k*T[ii]*N[ii]+np.array(results_E[ii]), yerr=results_uE[ii], fmt=".")
    # plt.title("Estymator energii calkowitej dla ukladu z potencjalem harmonicznym")
    # plt.xlabel("T")
    # plt.ylabel("<E>")
    # plt.legend(leg_entries)
    # plt.savefig("potencjal_harmoniczny.png")
    # plt.show()
    #calculations and plot for JL potential
    results_E = []
    results_uE = []
    for n in N:
        E = []
        uE = []
        for t in T:
            r = generate_system(n, 0., 0., 0., 1.)
            e, ue = dynamics(r, steps, sweep, sigma0, k, t, m, w, True, lj_eps, lj_sigma)
            E.append(e)
            uE.append(ue)
        results_E.append(E)
        results_uE.append(uE)
        steps = steps*3
    for ii in range(0, len(results_E)):
        plt.errorbar(T, 3. / 2. * k * T[ii] * N[ii] + np.array(results_E[ii]), yerr=results_uE[ii], fmt=".")
    plt.title("Estymator energii calkowitej dla ukladu z potencjalem Lennarda-Jonesa")
    plt.xlabel("T")
    plt.ylabel("<E>")
    plt.legend(leg_entries)
    plt.savefig("potencjal_LJ.png")

main()

"""
zmieniamy wzor na potencjal
potencjal Lennard-Jonesa (wzor jest na kartce z opisem doswiadczenia specjalistycznego)
sa dwie stale zalezne od rodzaju gazu, wybierzmy ze epsilon rowne 1 i rm=2


zrobic wykres E od T (punktowy z niepewnosciami E)
dla roznej liczby czastek (blad chcemy miec nie wiekszy niz 10% estymowanej wartosci)

1)
nieoddzialujacy gaz - pokazac ze dziala, wykres (E_calkowitej), omowic, (test chi^2?)
2)
rownanie dla lennarda Jonesa, wykres (E_calkowitej), omowic

"""