import os
import sys
sys.path.append(str(os.path.dirname(os.path.abspath(__file__)))+'\\..\\')
sys.path.append(str(os.path.dirname(os.path.abspath(__file__)))+'\\..\\..\\')
import createBeam as cb
from SAMPL.sourceCode.SAMPLcore.Particles import Electron
from SAMPL.sourceCode.SAMPLcore.SAMPLlab import Beamline
from SAMPL.sourceCode.SAMPLcore.SAMPLlab import Beam
from SAMPL.sourceCode.SAMPLcore.SAMPLlab import PhysicalUnits
import SAMPL.sourceCode.SAMPLcore.Components.Drift as D
import numpy as np


import matplotlib.pyplot as plt

beam1 = Beam.Beam(species=Electron.Electron, energy=4.0 * PhysicalUnits.MeV)
createBeam = cb.createBeam()
beam = createBeam.guassian(x=0.0, y=0.0, number=1000,
                           sigmaX=0.001, sigmaY=0.001)


a = D.Drift(name='', length=0.1)
n = 15
sX = np.zeros(n)
plt.scatter(beam.x, beam.y, s=10, c='r', alpha=0.5)
beam.bunchcharge=1000
for i in range(n):

    a.TrackSpaceCharge(beam)
    sX[i] = np.std(a.lastTrackedBeam.particles[:][0])

#print beam.particles


plt.scatter(beam.x, beam.y, s=10, c='b', alpha=0.5)

plt.show()

plt.plot(sX)

plt.show()
