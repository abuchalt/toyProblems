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
# ------------------------------------------------------------------------------

## Computational Parameters
# ------------------------------------------------------------------------------
phi = 1 # neutron flux [cm^-2.s^-1]

A1 = dm.Isotope(1, 1, 1, 0, 1)
A2 = dm.Isotope(2, 2, 0, 0, 1)
A3 = dm.Isotope(3, 3, 0, 0, 1)
A4 = dm.Isotope(4, 4, 0, 0, 1)
A5 = dm.Isotope(5, 5, 0, 0, 0)
A1.addDaughter(A2, 1)
A2.addDaughter(A3, 1)
A3.addDaughter(A4, 1)
A4.addDaughter(A5, 1)

isotopeList = [A1, A2, A3, A4, A5]

mySeries = dm.Series('mySeries', isotopeList)
# ------------------------------------------------------------------------------

# SIMULATION STUFF (WIP)
# ======================================================================

# Computation Parameters
# ----------------------------------------------------------------------
timeStop = 10 # [s]
deltaTime = timeStop/1e4 # [s]
timeDomain = np.arange(0, timeStop, deltaTime)

# Simulation Begin
# ----------------------------------------------------------------------

# Init
NA1 = [A1.N]
NA2 = [A2.N]
NA3 = [A3.N]
NA4 = [A4.N]
NA5 = [A5.N]

# Loop thru time
for i in range(len(timeDomain)-1):
    N_Δt = mySeries.timeStep(deltaTime, phi)
    NA1.append(N_Δt[0])
    NA2.append(N_Δt[1])
    NA3.append(N_Δt[2])
    NA4.append(N_Δt[3])
    NA5.append(N_Δt[4])

# Plotting
# ----------------------------------------------------------------------

# Isotopic Concentrations
plt.plot(timeDomain, NA1, label='NA1')
plt.plot(timeDomain, NA2, label='NA2')
plt.plot(timeDomain, NA3, label='NA3')
plt.plot(timeDomain, NA4, label='NA3')
plt.plot(timeDomain, NA5, label='NA3')
plt.xlabel('Time (s)')
plt.ylabel('Isotopic Abundance (mols/cc)')
plt.title('Bateman Problem')
plt.legend()
plt.savefig('results//IsotopicConcentrations.png')
plt.show()
plt.close()