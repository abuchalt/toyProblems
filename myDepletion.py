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
# ------------------------------------------------------------------------------

## Computational Parameters
# ------------------------------------------------------------------------------
phi = 1E14 # neutron flux [cm^-2.s^-1]

A1 = dm.Isotope(1, 1, 1, 0, 1)
A2 = dm.Isotope(1, 1, 0, 0, 1)
A3 = dm.Isotope(1, 1, 0, 0, 0)
A1.addDaughter(A2, 1)
A2.addDaughter(A3, 1)

isotopeList = [A1, A2, A3]

mySeries = dm.Series('mySeries', isotopeList)
# ------------------------------------------------------------------------------

# SIMULATION STUFF (WIP)
# ======================================================================

# Computation Parameters
# ----------------------------------------------------------------------
stopIrrad = 3*60*60*24 # Irradiate 3 days
timeStop = stopIrrad*2
deltaTime = timeStop/1e5
timeDomain = np.arange(0, timeStop, deltaTime)
timeDomainPlot = timeDomain/(60*60*24) # time in days
flux = 1e11
stopIndex = np.where(timeDomain>stopIrrad)[0][0]
fluxArray = np.zeros(len(timeDomain))
fluxArray[0:stopIndex] = flux

# Simulation Begin
# ----------------------------------------------------------------------

# Init
massH1 = []
massC12 = []
massCl35 = []
massCl37 = []
activityH3 = []
activityC14 = []
activityCl36 = []
activityCl38 = []
totalActivityArray = []

# Loop thru time
for i in range(len(timeDomain)):
    massH1.append(H1.N)
    massC12.append(C12.N)
    massCl35.append(Cl35.N)
    massCl37.append(Cl37.N)
    activityH3.append(H3.getActivity())
    activityC14.append(C14.getActivity())
    activityCl36.append(Cl36.getActivity())
    activityCl38.append(Cl38.getActivity())
    totalActivity = 0
    for item in isotopeList:
        totalActivity += item.getActivity()
        item.timeStep(deltaTime, fluxArray[i])
    totalActivityArray.append(totalActivity)

# Plotting
# ----------------------------------------------------------------------

# Isotopic Concentrations
plt.plot(timeDomainPlot, massH1, label='Hydrogen-1')
plt.plot(timeDomainPlot, massC12, label='Carbon-12')
plt.plot(timeDomainPlot, massCl35, label='Chlorine-35')
plt.plot(timeDomainPlot, massCl37, label='Chlorine-37')
plt.xlabel('Time (days)')
plt.ylabel('Isotopic Abundance (mols)')
plt.title('Irradiation of 1-mol Sample of PVC')
plt.legend()
plt.savefig('IsotopicConcentrations.png')
plt.show()
plt.close()