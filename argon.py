import numpy as np
from sys import argv
from time import time
from math import sqrt, pi
from numba import njit


def calculate_positions(a, n, N):
    atoms = np.zeros((N, 3))
    b0=np.array([a, 0, 0])
    b1=np.array([a/2, a/2*sqrt(3), 0])
    b2=np.array([a/2, a/6*sqrt(3), a*sqrt(2/3)])

    for i0 in range(n):
        for i1 in range(n):
            for i2 in range(n):
                i = i0 + i1*n + i2*n*n
                atoms[i] = (i0 - (n-1)/2)*b0 + (i1 - (n-1)/2)*b1 + (i2-(n-1)/2)*b2

    # np.savetxt('out_pos.txt', atoms, delimiter = '\t')
    return atoms

def calculate_momenta(k, T_0, N, m):
    energy = -1./2*k*T_0*np.log(np.random.uniform(size=(N, 3)))
    momenta = np.random.choice([-1., 1.], size=(N, 3))*np.sqrt(energy*2*m)
    momSum = np.sum(momenta, axis = 0)
    momenta -= momSum/N

    #np.savetxt('out_mom.txt', momenta, delimiter = '\t')
    return momenta

@njit #  = @jit(nopython=True)
def calculate_VFP(atoms, N, L, f, e, R):

    V = 0
    P = 0
    F_S = np.zeros((N, 3))
    F_P = np.zeros((N, N, 3))

    for i in range(N):
        ri = np.linalg.norm(atoms[i])
        if ri >= L:
            V += 1./2*f*(ri - L)**2
            F_S[i] = f*(L - ri)*atoms[i]/ri
        for j in range(i):
            rij = np.linalg.norm(atoms[i]-atoms[j])
            V += e*((R/rij)**12 - 2 * (R/rij)**6)
            F_P[i, j] = 12 * e * ((R/rij)**12 - (R/rij)**6) * ((atoms[i] - atoms[j])/(rij**2))
            F_P[j, i] = - F_P[i, j]

    F = np.sum(F_P, axis=1) + F_S
    P = np.sum(np.sqrt(np.sum(F_S**2, axis = 1)), axis=0)/(4*pi*L**2)

    # np.savetxt('out_for.txt', F, delimiter = '\t')
    return V, F, P


def simulation(atoms, momenta, N, L, f, e, R, tau, m, k, So, Sd, S_out, S_xyz):
    V, F, P = calculate_VFP(atoms, N, L, f, e, R)
    Tavr = Havr = Pavr = 0
    with open(argv[2], 'w+') as out, open(argv[3], 'w+') as outxyz:
        out.write("t \t\t H \t\t V \t\t T \t\t P \n")
        for s in range(1, So + Sd + 1):

                momenta = momenta + 1./2*F*tau
                atoms = atoms + 1./m * momenta*tau
                V, F, P = calculate_VFP(atoms, N, L, f, e, R)
                momenta = momenta + 1./2* F*tau

                Ekin = np.linalg.norm(momenta, axis = 1)**2 / (2*m)
                Ekin_sum = np.sum(Ekin, axis = 0)
                H = Ekin_sum + V
                T= 2 / (3*N*k)*Ekin_sum

                if(s % S_out == 0): #TODO
                    # zapis chwilowych charakterystyk t, H, V, T, P w formacie wykładniczym do out.txt
                    out.write("{:.3f}".format(s*tau) + "\t")
                    out.write("{:.5E}".format(H) + "\t")
                    out.write("{:.5E}".format(V) + "\t")
                    out.write("{:.5E}".format(T) + "\t")
                    out.write("{:.5E}".format(P) + "\n")

                if(s % S_xyz == 0):
                    # zapis wsp. atomów oraz ich Ekin do avs.dat: xi, yi, zi, Ekini
                    data_xyz = np.concatenate((atoms, Ekin[:, None]), axis = 1)
                    np.savetxt(outxyz, data_xyz, delimiter = '\t')
                    outxyz.write('\n\n')

                if(s >= So):
                    Tavr += T
                    Pavr += P
                    Havr += H

    Tavr /= Sd
    Pavr /= Sd
    Havr /= Sd
    print("Tavr: ", Tavr)
    print("Pavr: ", Pavr)
    print("Havr: ", Havr)


if __name__ == "__main__":

    tic = time()

    params = {}
    with open(argv[1]) as f:
        for line in f:
            val, key = line.split("#")
            params[key.strip()] = float(val.strip())
    params['N'] = int(params['n']*params['n']*params['n'])

    atoms = calculate_positions(params['a'], int(params['n']), params['N'])
    momenta = calculate_momenta(params['k'], params['T_0'], params['N'], params['m'])
    # V, F, P = calculate_VFP(atoms, params['N'], params['L'], params['f'], params['e'], params['R'])
    simulation(atoms, momenta, params['N'], params['L'], params['f'], params['e'], params['R'], params['tau'], params['m'], params['k'], int(params['So']), int(params['Sd']), int(params['S_out']), int(params['S_xyz']))

    toc = time()
    print('execution time: ', (toc-tic))
