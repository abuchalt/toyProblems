## Depletion Code
# ------------------------------------------------------------------------------
# Adam Buchalter
#
# Machine-accurate evaluation of a decay series with production terms
# ------------------------------------------------------------------------------

## Imports
# ------------------------------------------------------------------------------
import numpy as np
from scipy import linalg
import depletionModule as dm
import matplotlib.pyplot as plt
import pandas as pd
# ------------------------------------------------------------------------------

# SIMULATION STUFF
# ======================================================================

# Computation Parameters
# ----------------------------------------------------------------------
radStop = 3*3.154e7 # [a->s]
timeStop = 30*3.154e7 # [a->s]
deltaTime = timeStop/1e2 # [s]
timeDomain = np.arange(0, timeStop, deltaTime)
phi = 1E13 # neutron flux [cm^-2.s^-1]

# Uranium Fuel
# ----------------------------------------------------------------------

# Init Isotopes
Enrichment = 0.05
M_UO2 = 270.03 # [g/mol]
ρ_UO2 = 10.97 # [g/cc]
N_UO2 = ρ_UO2/M_UO2 # [mol/cc]
sigma_nyU235 = 0.08741 # [barns] JAEA https://wwwndc.jaea.go.jp/cgi-bin/Tab80WWW.cgi?lib=J40&iso=U235
sigma_fU235 = 1.218 # [barns] JAEA
λ_U235 = np.log(2)/(7.040E8*3.154E7) # [s^-1]
sigma_nyU238 = 0.07016 # [barns] JAEA
sigma_fU238 = 0.3064 # [barns] JAEA
λ_U238 = np.log(2)/(4.463E9*3.154E7) # [s^-1]

U235 = dm.Isotope(235, 92, N_UO2*Enrichment, sigma_nyU235, λ_U235, sigma_fU235)
U238 = dm.Isotope(238, 92, N_UO2*(1-Enrichment), sigma_nyU238, λ_U238, sigma_fU238)
isotopeList = [U235,U238]

prodsData = pd.read_csv('data\\fisprods.csv')
NProds = len(prodsData)
for i in range(NProds):
    thisProd = prodsData.iloc[i]
    Z = thisProd['Z']
    A = thisProd['A']
    thisfrac235 = thisProd['235U']
    thisfrac238 = thisProd['238U']
    thislambda = thisProd['lambda']
    if thisfrac235 > 0 or thisfrac238 > 0:
        thisIsotope = dm.Isotope(A, Z, 0, 0, thislambda, 0)
        if thisfrac235 > 0:
            U235.addFProduct(thisIsotope, thisfrac235)
        if thisfrac238 > 0:
            U238.addFProduct(thisIsotope, thisfrac238)
        isotopeList.append(thisIsotope)

mySeries = dm.Series('mySeries', isotopeList)
NIsotopes = len(isotopeList)
NTimeSteps = len(timeDomain)

# Simulation Begin
# ----------------------------------------------------------------------

# Init
outData = np.zeros(shape=(NIsotopes,NTimeSteps))
for i in range(NIsotopes):
    outData[i,0] = isotopeList[i].N

# Loop thru time
for i in range(len(timeDomain)-1):
    if timeDomain[i] < radStop:
        N_Δt = mySeries.timeStep(deltaTime, phi)
    else:
        N_Δt = mySeries.timeStep(deltaTime, 0)
    for j in range(NIsotopes):
        outData[j,i+1] = N_Δt[j]


# Plotting
# ----------------------------------------------------------------------

# Isotopic Concentrations
for i in range(3, NIsotopes):
    myLabel = isotopeList[i].getZaid()
    plt.plot(timeDomain/3.154e7, outData[i], label=myLabel)
plt.xlabel('Time (a)')
plt.ylabel('Isotopic Abundance (mols/cc)')
plt.title('Depletion in UO2 Fuel')
plt.legend()
plt.savefig('results//DepletionFig.png')
plt.show()
plt.close()

## Computational Parameters
# ------------------------------------------------------------------------------
# phi = 1 # neutron flux [cm^-2.s^-1]

# A1 = dm.Isotope(1, 1, 1, 0, 1)
# A2 = dm.Isotope(2, 2, 0, 0, 1)
# A3 = dm.Isotope(3, 3, 0, 0, 1)
# A4 = dm.Isotope(4, 4, 0, 0, 1)
# A5 = dm.Isotope(5, 5, 0, 0, 0)
# A1.addDaughter(A2, 1)
# A2.addDaughter(A3, 1)
# A3.addDaughter(A4, 1)
# A4.addDaughter(A5, 1)

# isotopeList = [A1, A2, A3, A4, A5]

# mySeries = dm.Series('mySeries', isotopeList)
# ------------------------------------------------------------------------------

# Test Sim
# ----------------------------------------------------------------------
# # Init
# NA1 = [A1.N]
# NA2 = [A2.N]
# NA3 = [A3.N]
# NA4 = [A4.N]
# NA5 = [A5.N]

# # Loop thru time
# for i in range(len(timeDomain)-1):
#     N_Δt = mySeries.timeStep(deltaTime, phi)
#     NA1.append(N_Δt[0])
#     NA2.append(N_Δt[1])
#     NA3.append(N_Δt[2])
#     NA4.append(N_Δt[3])
#     NA5.append(N_Δt[4])

# # Plotting
# # ----------------------------------------------------------------------

# # Isotopic Concentrations
# plt.plot(timeDomain, NA1, label='NA1')
# plt.plot(timeDomain, NA2, label='NA2')
# plt.plot(timeDomain, NA3, label='NA3')
# plt.plot(timeDomain, NA4, label='NA3')
# plt.plot(timeDomain, NA5, label='NA3')
# plt.xlabel('Time (s)')
# plt.ylabel('Isotopic Abundance (mols/cc)')
# plt.title('Bateman Problem')
# plt.legend()
# plt.savefig('results//IsotopicConcentrations.png')
# plt.show()
# plt.close()