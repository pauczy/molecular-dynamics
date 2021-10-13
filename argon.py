import numpy as np
import sys
import math
from math import sqrt

class Parameters:
    def __init__(self):
        infile = sys.argv[1]
        lines = np.loadtxt(infile, comments = "#")
        self.n = int(lines[0])
        self.m = lines[1]
        self.e = lines[2]
        self.R = lines[3]
        self.f = lines[4]
        self.L = lines[5]
        self.a = lines[6]
        self.T_0 = lines[7]
        self.tau = lines[8]
        self.S_o = lines[9]
        self.S_d = lines[10]
        self.S_out = lines[11]
        self.S_xyz = lines[12]

        self.N = self.n*self.n*self.n
        self.k = 0.00831

def calculate_positions(atoms, parameters):
    b0=np.array([parameters.a, 0, 0])
    b1=np.array([parameters.a/2, parameters.a/2*sqrt(3), 0])
    b2=np.array([parameters.a/2, parameters.a/6*sqrt(3), parameters.a*sqrt(2/3)])

    for i0 in range(parameters.n):
        for i1 in range(parameters.n):
            for i2 in range(parameters.n):
                i = i0 + i1*parameters.n + i2*parameters.n*parameters.n
                atoms[i] = (i0 - (parameters.n-1)/2)*b0 + (i1 - (parameters.n-1)/2)*b1 + (i2-(parameters.n-1)/2)*b2

    np.savetxt('out_pos.txt', atoms, delimiter = '\t')

def calculate_momenta(parameters):
    energy = -1./2*parameters.k*parameters.T_0*np.log(np.random.uniform(size=(parameters.N, 3)))
    momenta = np.random.choice([-1., 1.], size=(parameters.N, 3))*np.sqrt(energy*2*parameters.m)
    momSum = np.sum(momenta, axis = 0)
    momenta -= momSum/parameters.N

    np.savetxt('out_mom.txt', momenta, delimiter = '\t')

#potencjal ogolnia wartosc, sily dla wszystkich atomow wypisac
def calculate_potential(atoms, parameters):
    V = 0.
    F = np.zeros((parameters.N, 3))
    P = 0

    #potencjal od scianek
    V_s = np.zeros((parameters.N))
    r = np.sqrt(np.sum(atoms**2, axis = 1))
    V_s[r >= parameters.L] = 1./2*parameters.f*np.square(r-parameters.L)[r >= parameters.L]
    V += np.sum(V_s, axis = 0)
    print(V_s)
    print('V_s ',  np.sum(V_s, axis=0))




    #sily odpychania od scianek
    F_s = np.zeros((parameters.N, 3))
    F_s_calc = np.repeat(parameters.f*(parameters.L-r)[np.newaxis, :], 3, 0).T * np.divide(atoms, np.repeat(r[np.newaxis, :], 3, 0).T)
    F_s[r >= parameters.L, :] = F_s_calc[r >= parameters.L, :]
    F += F_s

    #cisnienie chwilowe
    Fnorm_sum = np.sum(np.sum(F_s**2, axis = 1), axis=0)
    P = 1./(4*math.pi*parameters.L**2)*Fnorm_sum


    #potencjał  par atomowych
    V_vdw = np.zeros((parameters.N, parameters.N))
    rij = np.zeros((parameters.N, parameters.N))
    for i in range(parameters.N):
        for j in range(i):
                rij[i, j] = math.sqrt(np.sum(np.square(atoms[i] - atoms[j])))
                V_vdw[i, j] = parameters.e * ( (parameters.R/ rij[i,j])**12 - 2*(parameters.R/rij[i,j])**6 )
    #Rr = np.divide(parameters.R, rij)
    #V_vdw = parameters.e*(np.power(Rr, 12) - 2*np.power(Rr, 6))
    V += np.nansum(V_vdw)

    #sily par atomowych
    F_vdw = np.zeros((parameters.N, parameters.N, 3))
    for i in range(parameters.N):
        for j in range(i):
            rij2 = math.sqrt(np.sum(np.square(atoms[i] - atoms[j])))
            F_vdw[i, j, :] = 12*parameters.e*((parameters.R/rij2)**12 - (parameters.R/rij2)**6)*(atoms[i]-atoms[j])/rij2**2

            F[i, :] += F_vdw[i, j, :]
            F[j, :] -= F_vdw[i, j, :]

    print(V)

    np.savetxt('out_for.txt', F, delimiter = '\t')







if __name__ == "__main__":

    params = Parameters()
    atoms = np.zeros((params.N, 3))
    calculate_positions(atoms, params) #polozenia zapisane w atoms
    calculate_momenta(params)
    calculate_potential(atoms, params)
