import os
import sys
sys.path.append(str(os.path.dirname(os.path.abspath(__file__)))+'\\..\\')
sys.path.append(str(os.path.dirname(os.path.abspath(__file__)))+'\\..\\..\\')
from SAMPL.sourceCode.SAMPLcore.Components import Drift as d
from SAMPL.sourceCode.SAMPLcore.Components import Dipole as D
from SAMPL.sourceCode.SAMPLcore.Components import Quadrupole as Q
from SAMPL.sourceCode.SAMPLcore.Components import Screen as S
from SAMPL.sourceCode.SAMPLcore.Components import OrbitCorrector as C
from SAMPL.sourceCode.SAMPLcore.Components import BeamPositionMonitor as BPM
from SAMPL.sourceCode.SAMPLcore.SAMPLlab import Beamline
from SAMPL.sourceCode.SAMPLcore.SAMPLlab import Beam
from SAMPL.sourceCode.SAMPLcore.Particles import Electron
from SAMPL.sourceCode.SAMPLcore.SAMPLlab import PhysicalUnits
import SAMPL.sourceCode.SAMPLcore.Physics.ComputeTransferMatrix as CTM
import numpy as np


beam1 = Beam.Beam(species=Electron.Electron, energy=4.5 * PhysicalUnits.MeV)
ptcle1 = [0.001, 0, 0, 0, 0, 0]
ptcle2 = [0, 0, 0.001, 0, 0, 0]
ptcle3 = [0, 0, 0, 0, 0, 0]
beam1.particles = np.array([ptcle1, ptcle2, ptcle3])

Line = Beamline.Beamline(componentlist=[])
Line.componentlist.append(d.Drift(name='drift1', length=0.1))
Line.componentlist.append(Q.Quadrupole(name='quad1', gradient=0.0, length=0.05))

beam2 = Line.TrackMatlab([0, len(Line.componentlist) - 1], beam1)

traj = np.array([0, 0, 0, 0, 0, 0])
a, b = CTM.ComputeTransferMatrix(Line, [0, 1], beam1, traj)

print b
print a
