import numpy as np
import sys
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


if __name__ == "__main__":

    params = Parameters()
    atoms = np.zeros((params.N, 3))
    calculate_positions(atoms, params)
    calculate_momenta(params)
