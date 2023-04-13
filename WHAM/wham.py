import numpy as np
import sys 
import numpy as np 
import pytraj as pt
import matplotlib.pyplot as plt
import pandas as pd

from string import *
from math import *
from array import *

#---------------------------------

Nwindows = 12

#-----------------calc UBias:
def biasPotential(dihedral, barrier1, phase1, barrier2, phase2):
    ''' '''
    dihedral = np.radians(dihedral)

    original = 4 * (2.5 *(1 + np.cos(2 * dihedral - np.pi)))
    
    return ((barrier1 *(1 + np.cos(dihedral - phase1))) +  4 * (barrier2 *(1 + np.cos(2 * dihedral - phase2)))) - original


barrierHeightTwoFold = [0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 15.0, 18.0, 20.0]
barrierParmedTwoFold = [float(a / 8.0) for a in barrierHeightTwoFold]
Barrier1 = [1 for rep in range(Nwindows)]
Phase1 = [np.pi for rep in range(Nwindows)]
Barrier2 = barrierParmedTwoFold
Phase2 = [np.pi for rep in range(Nwindows)]


def WHAMresults(histogramlabel, outputResults):

    Nwindows = 12 # number of windows
    windows = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] # indices of the replicas
    Nbins = 120 # number of bins in the input histograms (i.e. each data file has Nbins entries)
    Ni = 0.0 # initialize count number to 0 (see below)
    kBT = 0.001987204 * 300 # in kcal/mol

    MAXITER = 20000 # max number of iterations
    counts = []
    UBiasPotential = []
    # arrays for the pmf
    bins = []
    prob = []
    F = []
    PMF = []

    #initiate arrays:

    for i in range(Nwindows):	#initial guess for F is zero
        F.append(0.0)
        counts.append([])
        UBiasPotential.append([])
        for j in range(Nbins):	
            counts[i].append(0)	#set initials counts to zero
            UBiasPotential[i].append(0)
            bins.append(0) #too much?
            prob.append(0)
            PMF.append(0)

    #------------readin datFiles
    for i in range(Nwindows):
        stringvalue = f'{windows[i]}'
        #print(stringvalue)
        # read in the datafiles -- data0.dat, data1.dat ....
        with open(histogramlabel + stringvalue + ".dat", "r")	as f:
            for j in range(Nbins):
                line = f.readline()
                s = line.split()
                #counts[i][j] = 0
                try:
                  counts[i][j] = float(s[1])
                except:
                    print(j, s[1])
                if i == 0:
                  bins[j] = float(s[0])	#read bins
                  Ni = Ni + counts[i][j]	#calculate number of counts per simulation; these were read from file above

    for i in range(Nwindows):	#calculate UBiasPotential for all simulations and bins
        for j in range(Nbins): # j loop dihedral angles
            UBiasPotential[i][j] = biasPotential(bins[j], Barrier1[i], Phase1[i], Barrier2[i], Phase2[i])


    for k in range(MAXITER):	#repeat calculation of prob and F for MAXITER
        for j in range(Nbins):
            denominator = 0.0
            numerator = 0.0
            for i in range(Nwindows):	#calculation of P(x)
                denominator = denominator + (Ni * exp((F[i] - UBiasPotential[i][j]) / kBT))
                numerator = numerator + counts[i][j]
            prob[j] =  numerator / denominator

        for i in range(Nwindows):	#calculation of F[i]
            fprob = 0
            for j in range(Nbins):
                fprob = fprob + prob[j] * exp(-UBiasPotential[i][j] / kBT)
            F[i] = -kBT * log(fprob)

    with open(outputResults, 'w') as g:
        for j in range(Nbins):
            if prob[j] == 0:
                PMF[j] = 0
            else:
                PMF[j] = -kBT *log(prob[j])
            g.write(str(bins[j]) + "     " + str(PMF[j]) + "\n")


WHAMresults('histogram', 'result')