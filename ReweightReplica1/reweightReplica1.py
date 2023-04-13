import numpy as np
import pickle
from scipy.special import i0
import mpmath


def vonmisesKDE(data, kappa):
       ''' 
       get vonMises KDE on xgrid
       '''
       
       kde = np.exp(kappa*np.cos(xgrid[:, None]- data[None, :])).sum(1)/np.float128((2*np.pi*mpmath.besseli(0,kappa)))
       kde /= len(data) # Normalization
       return kde


def biasPotential(dihedral, barrier1, phase1):
       '''
       values of the oneFold bias potential in replica1
       '''
       
       return barrier1 *(1 + np.cos(dihedral - phase1))


def reweight(Pbiased, bias):
       '''
       reweighting of biased probabilities 
       '''

       return np.exp(bias/(k_b*T)) * Pbiased


def freeE(Probability):
       '''
       Free Energy calculation according to unbiased probabilities
       '''
       
       return -np.log(Probability)* k_b * T

#### MAIN ####

k_b = 0.001987204 # k in kcal/mol*K
T = 300 # temperature in K
n_gridpoints = 100
xgrid = np.linspace(-np.pi, np.pi, n_gridpoints) # grid for kernel density estimation


omegaValues = pickle.load(open('exampleOmega.p', 'rb')) #dummy data for sampling in replica1

KDE = vonmisesKDE(omegaValues, 300)
reweightedKDE = [reweight(kde, bias) for kde, bias in  zip(KDE, biasPotential(xgrid, barrier1=1, phase1=np.pi))] #exemplary bias in replica1 
reweightedKDE /= np.trapz(reweightedKDE, xgrid) #Normalization

G = freeE(reweightedKDE)

