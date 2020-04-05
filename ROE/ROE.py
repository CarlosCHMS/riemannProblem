# Programed in Python3 but run in python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 15:39:49 2020

@author: carlos
"""


import numpy as np
import matplotlib.pyplot as plt
import subprocess as sp
import characteristics as ch


def trapz(x, y):

    # Integration by trapezoidal method
    ans = 0
    for ii in range(0, len(x)-1):
        ans += (x[ii+1] - x[ii])*(y[ii] + y[ii+1])/2

    return ans


class configuration():

    def __init__(self, Np=100, xtot=50, dt=0.0001, Ns=100,
                 xchange=25, p0=1e6, T0=300, p1=1e5, T1=300, nSave=4, R=287,
                 gamma=1.4, flux='ROE'):

        self.var = dict()

        self.var['Np'] = Np
        self.var['xtot'] = xtot
        self.var['dt'] = dt
        self.var['Ns'] = Ns

        self.var['xchange'] = xchange
        self.var['p0'] = p0
        self.var['T0'] = T0
        self.var['p1'] = p1
        self.var['T1'] = T1

        self.var['nSave'] = nSave

        self.var['R'] = R
        self.var['gamma'] = gamma

        self.var['flux'] = flux

        # Calculated values:
        self.var['dx'] = xtot/(Np-1)
        self.var['n'] = xchange/self.var['dx']
        self.var['cc'] = self.var['dt']/self.var['dx']
        self.var['saveStep'] = Ns/nSave

        self.var['fluxType'] = 0
        if self.var['flux'] == 'ROE':
            self.var['fluxType'] = 1
        elif self.var['flux'] == 'AUSM':
            self.var['fluxType'] = 2

        self.x = xtot*np.arange(0, Np)/Np

    def printVar(self):

        for n in self.var.keys():
            print(n, ' = ', self.var[n])

        return None

    def printVarAuxFile(self):

        for v in ['Ns', 'Np']:
            print('%s, %i,' % (v, self.var[v]))

        for v in ['n', 'saveStep', 'p0', 'T0', 'p1', 'T1', 'cc', 'R', 'gamma']:
            print('%s, %f,' % (v, self.var[v]))

        for v in ['fluxType']:
            print('%s, %i,' % (v, self.var[v]))

        return None

    def writeAuxFile(self):

        f = open('aux.csv', 'w')

        for v in ['Ns', 'Np']:
            f.write('%s, %i,\n' % (v, self.var[v]))

        for v in ['n', 'saveStep', 'p0', 'T0', 'p1', 'T1', 'cc', 'R', 'gamma']:
            f.write('%s, %f,\n' % (v, self.var[v]))

        for v in ['fluxType']:
            f.write('%s, %i,\n' % (v, self.var[v]))

        f.close()

        return None


class output():

    def __init__(self, conf):

        self.output = open("./output.csv", 'r')

        self.R = conf.var['R']
        self.g = conf.var['gamma']

        self.x = conf.x
        self.dt = conf.var['dt']
        self.Np = conf.var['Np']

        self.iSol = []
        self.CFL = []
        self.sol = []
        self.Ns = 0

        self.readOutput()

        self.varList = ['t', 'x', 'rho', 'm', 'e', 'u', 'p', 'T']
        self.varUnit = {'t': '[s]', 'x': '[m]', 'rho': '[kg/m3]',
                        'm': '[kg/m2 s]', 'e': '[J/m3]', 'u': '[m/s]',
                        'p': '[Pa]', 'T': '[K]'}

    def readOutput(self):

        ii = 0
        for row in self.output:
            aux = row.split(',')

            if ii == 0:
                # print(ii, kk, aux)
                N = int(aux[0])
                self.iSol.append(int(aux[1]))
                self.CFL.append(float(aux[2]))
                U = np.zeros((N, 3))
            else:
                for jj in range(0, N):
                    U[jj, ii-1] = float(aux[jj])

            ii += 1
            if ii == 4:
                ii = 0
                self.sol.append(U.copy())

        self.Ns = len(self.sol)

        return None

    def calcVar(self, n):

        # Initialize variables
        var = dict()
        for v in self.varList:
            if v not in ['t', 'x']:
                var[v] = np.zeros((self.Np,))

        var['t'] = self.iSol[n]*self.dt
        var['x'] = self.x

        for jj in range(0, self.Np):

            rho = self.sol[n][jj, 0]
            m = self.sol[n][jj, 1]
            e = self.sol[n][jj, 2]
            u = m/rho
            p = (e - m*u/2)*(self.g-1)
            T = p/(self.R*rho)

            var['rho'][jj] = rho
            var['m'][jj] = m
            var['e'][jj] = e

            var['u'][jj] = u
            var['p'][jj] = p
            var['T'][jj] = T

        var['unit'] = self.varUnit

        return var

    def integrals(self):

        for ii in range(0, self.Ns):

            t = self.iSol[ii]*self.dt
            integ = trapz(self.x, self.sol[ii][:, 0])
            print("\n rho integral for t = %.3f: %.2f" % (t, integ))

        for ii in range(0, self.Ns):

            t = self.iSol[ii]*self.dt
            integ = trapz(self.x, self.sol[ii][:, 2])
            print("\n e integral for t = %.3f: %.0f" % (t, integ))

        return None


if __name__ == "__main__":

    conf1 = configuration(Ns=300, Np=500, xtot=50,
                          xchange=25, p0=1e6, T0=300, p1=1e5, T1=300,
                          dt=0.00013, flux='ROE')

    print("\nSolve Riemann problem with a CFD solver:")
    conf1.printVar()
    print("\nAux File:")
    conf1.printVarAuxFile()
    conf1.writeAuxFile()
    print()

    sp.call("./bin/Debug/ROE", shell=True)

    out1 = output(conf1)
    print("CFL: ", out1.CFL)

    cfdVar = out1.calcVar(out1.Ns-1)

    # Solving using characteristics
    char1 = ch.problem(xchange=25.0, p1=1e5, T1=300, p4=1e6, T4=300)
    charVar = char1.calcVar(cfdVar['t'], out1.x)

    plt.close('all')

    v = 'rho'
    plt.figure()

    # ploting cfd results
    plt.plot(cfdVar['x'], cfdVar[v])

    # ploting characteristics results
    plt.plot(charVar['x'], charVar[v], '--')

    plt.legend(['CFD', 'charac.'])

    plt.xlabel("%s %s" % ('x', cfdVar['unit']['x']))
    plt.ylabel("%s %s" % (v, cfdVar['unit'][v]))

    title = "t: %f [s]" % cfdVar['t']
    plt.title(title)

    plt.grid(True)
    plt.show()

    out1.integrals()
    #print(trapz(charVar['x'], charVar['e']))
