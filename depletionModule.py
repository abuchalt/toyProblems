## Depletion Module
# ------------------------------------------------------------------------------
# Adam Buchalter
#
# Module for Building Depletion Series
# ------------------------------------------------------------------------------

## Imports
# ------------------------------------------------------------------------------
import numpy as np
from scipy import linalg
# ------------------------------------------------------------------------------

## Global Vars
# ------------------------------------------------------------------------------
N_A = 6.022e23 # Avogadro's Number (mol^-1)
unit = 1E-24 # Unit Conversion (cm^2/barn)
# ------------------------------------------------------------------------------

## Classes
# ------------------------------------------------------------------------------
class Isotope:
    def __init__(self, A:int, Z:int, N:float=0.0, sigma_ny=0.0, decayConst=0.0, sigma_f=0.0):
        self.A = A # Mass Number
        self.Z = Z # Atomic Number
        self.N = N # Number density of isotope [mols/cm^-3]
        self.zaid = Z*1000 + A # zaid (ZZZAAA)
        self.sigma_ny = sigma_ny # transmutation (n,gamma) cross-section [barns]
        self.sigma_f = sigma_f
        self.decayConst = decayConst # decay constant [s^-1]
        self.daughters = [] # Decay Daughters
        self.products = [] # n,gamma Products
        self.fproducts = [] # fission Products

    def __str__(self):
        return f"{self.N}"
    
    def getZaid(self):
        return f"{self.zaid}"
    
    def addDaughter(self, daughter:"Isotope", fraction:float):
        self.daughters.append((daughter, fraction))
    
    def addProduct(self, product:"Isotope", fraction:float):
        self.products.append((product, fraction))

    def addFProduct(self, product:"Isotope", fraction:float):
        self.fproducts.append((product, fraction))

    def timeStep(self, deltaT=1, flux=0):
        # Neutron Absorption
        deltaFluence = flux*deltaT # Integrate flux to produce fluence
        reactionN = self.N*self.sigma_ny*deltaFluence # # of reactions in mols
        self.N -= reactionN # Remove from mols of isotope
        for tuple in self.products: # Add to mols of products
            tuple[0].N += reactionN*tuple[1]
        # Decay
        decayN = self.N*self.decayConst*deltaT # # of decays in mols
        self.N -= decayN # Remove from mols of isotope
        for tuple in self.daughters: # Add to mols of products
            tuple[0].N += decayN*tuple[1]

    def getActivity(self):
        #decayRate = self.N*self.decayConst*N_A # Decay Rate in becquerel
        #decayRate = self.N*self.decayConst*N_A/37000000000 # Decay Rate in curie
        decayRate = self.N*self.decayConst*N_A/37000 # Decay Rate in microcurie
        return decayRate
    
class Series:
    def __init__(self, name:str, isotopes=[]):
        self.name = name # Name
        self.isotopes = isotopes # List of Isotopes in Series

    def buildMat(self, phi):
        n = len(self.isotopes) # number of isotopes in series
        A = np.zeros((n, n), dtype=float)
        N_0 = np.zeros(n, dtype=float)
        i = 0
        for isotope in self.isotopes:
            thisLambda = isotope.decayConst
            thisSigma = isotope.sigma_ny
            thisFSigma = isotope.sigma_f
            thisN = isotope.N
            N_0[i] = thisN
            A[i, i] += - thisLambda - phi*(thisSigma + thisFSigma)*unit
            for daughter in isotope.daughters:
                # daughterzaid = daughter[0].zaid
                fraction = daughter[1]
                j = self.isotopes.index(daughter[0])
                A[j, i] += fraction*thisLambda
            for product in isotope.products:
                # daughterzaid = daughter[0].zaid
                fraction = product[1]
                j = self.isotopes.index(product[0])
                A[j, i] += fraction*phi*thisSigma*unit
            for fproduct in isotope.fproducts:
                fraction = fproduct[1]
                j = self.isotopes.index(fproduct[0])
                A[j, i] += fraction*phi*thisFSigma*unit
            i += 1
        return A, N_0

    def timeStep(self, Δt, phi):
        A, N_0 = self.buildMat(phi)
        At = A*Δt
        N_Δt = linalg.expm(At) @ N_0
        i = 0
        for isotope in self.isotopes:
            isotope.N = N_Δt[i]
            i+=1
        return N_Δt


# ------------------------------------------------------------------------------