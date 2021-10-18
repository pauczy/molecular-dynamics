import numpy as np
import sys
import time
import math
from math import sqrt


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

    np.savetxt('out_pos.txt', atoms, delimiter = '\t')
    return atoms

def calculate_momenta(k, T_0, N, m):
    energy = -1./2*k*T_0*np.log(np.random.uniform(size=(N, 3)))
    momenta = np.random.choice([-1., 1.], size=(N, 3))*np.sqrt(energy*2*m)
    momSum = np.sum(momenta, axis = 0)
    momenta -= momSum/N

    np.savetxt('out_mom.txt', momenta, delimiter = '\t')
    return momenta

#potencjal ogolnia wartosc, sily dla wszystkich atomow wypisac
def calculate_VFP(atoms, N, L, f, e, R):
    V = 0.
    F = np.zeros((N, 3))

    #wall potential
    V_s = np.zeros((N))
    r = np.sqrt(np.sum(atoms**2, axis = 1))
    V_s[r >= L] = 1./2*f*np.square(r-L)[r >= L]
    V += np.sum(V_s, axis = 0)

    #wall forces
    F_s = np.zeros((N, 3))
    F_s[r >= L] =f*(L-r[r >= L, np.newaxis])*atoms[r >= L]/r[r >= L, np.newaxis]
    F += F_s

    #van der waals potential + forces between atom pairs
    V_vdw = np.zeros((N, N))
    F_vdw = np.zeros((N, N, 3))
    for i in range(N):
        for j in range(i):
            rij = np.linalg.norm(atoms[i] - atoms[j])
            V_vdw[i, j] = e * ( (R/ rij)**12 - 2*(R/rij)**6 )
            F_vdw[i, j, :] = 12*e*((R/rij)**12 - (R/rij)**6)*(atoms[i]-atoms[j])/rij**2
            F[i, :] += F_vdw[i, j, :]
            F[j, :] -= F_vdw[i, j, :]
    V += np.sum(V_vdw)

    #temporary pressure
    P = np.sum(np.sum(F_s**2, axis = 1), axis=0)/(4*math.pi*L**2)

    print('V_s: ',  np.sum(V_s, axis=0))
    print('V: ', V)
    print('P: ', P)
    np.savetxt('out_for.txt', F, delimiter = '\t')

    return V, F, P

def simulation(atoms, momenta, N, L, f, e, R, tau, m, k, So, Sd):
    V, F, P = calculate_VFP(atoms, N, L, f, e, R)
    Tavr = Havr = Pavr = 0
    with open(sys.argv[2], 'w+') as outfile, open(sys.argv[3], 'w+') as outxyz:
        for s in range(1, So + Sd + 1):
                momenta = momenta + 1./2*F*tau
                atoms = atoms + 1./m * momenta*tau
                V, F, P = calculate_VFP(atoms, N, L, f, e, R)
                momenta = momenta + 1./2* F*tau

                Ekin = np.sum(np.linalg.norm(momenta)**2 / (2*m), axis=0)
                H = Ekin + V
                T= 2 / (3*N*k)*Ekin

                if(s % S_out == 0):

                if(s % S_xyz == 0):

                if(s >= So):





if __name__ == "__main__":

    tic = time.time()

    params = {}
    with open(sys.argv[1]) as f:
        for line in f:
            val, key = line.split("#")
            params[key.strip()] = float(val.strip())
    params['N'] = int(params['n']*params['n']*params['n'])
    params['k'] = 0.00831

    atoms = calculate_positions(params['a'], int(params['n']), params['N'])
    momenta = calculate_momenta(params['k'], params['T_0'], params['N'], params['m'])
    V, F, P = calculate_VFP(atoms, params['N'], params['L'], params['f'], params['e'], params['R'])
    simulation(atoms, momenta, params['N'], params['L'], params['f'], params['e'], params['R'], params['tau'], params['m'], params['k'], int(params['So']), int(params['Sd']))
    print(params['tau'])
    toc = time.time()
    print('execution time in ms: ', 1000*(toc-tic))
