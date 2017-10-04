# from SAMMcore.Components import *
from SAMMcore.Components import Drift
from SAMMcore.Components import Dipole
from SAMMcore.Components import Quadrupole
from SAMMcore.SAMMlab import Beamline
from SAMMcore.SAMMlab import Beam
from SAMMcore.Particles import Positron
from SAMMcore.SAMMlab import PhysicalUnits
import numpy as np

[t1, t2] = np.array([[1], [2]])

beam1 = Beam.Beam(species=Positron.Positron)

beam1.energy = 1.20 * PhysicalUnits.GeV

# help(comp)

drift1 = Drift.Drift(length=1.5)
quad1 = Quadrupole.Quadrupole(length=0.3, gradient=0.8)
dip1 = Dipole.Dipole(field=1, length=0.5, gradient=0.8)


dip1.theta = 0.2
dip1.curvature = 3
print('dip1.curvature = ', dip1.curvature)

print('dip1.curvature = ', dip1.curvature)

# print dip1.curvature

beamline1 = Beamline.Beamline(componentlist=[quad1, drift1])

# print beamline1.precision
# print beamline1.ComputePositions()

beam1 = Beam.Beam(species=Positron.Positron, energy=1.20 * PhysicalUnits.GeV)
# beam1 = Beam.Beam(species=Positron, energy = 2)


# print('beam1.energy = ', beam1.energy)
# print('beam1.rigidity = ',beam1.rigidity )
# print('beam1.momentum = ',beam1.momentum )

ptcle1 = [0.001, 0, 0, 0, 0, 0]
ptcle2 = [0, 0, 0.001, 0, 0, 0]
parts = np.array([ptcle1, ptcle2])

# print parts[:,0]
# print parts[:,1]
# print parts[:,2]
# print parts[:,3]

beam1.particles = np.array([ptcle1, ptcle2])
beam2 = beamline1.TrackMatlab([0, 2], beam1)
print("beam2 = ", beam2)

print np.array([[9.01351346e-04, -5.97789022e-05, 0.00000000e+00,
                0.00000000e+00, -2.50103285e-09, 0.00000000e+00],
                [0.00000000e+00, 0.00000000e+00, 1.09921487e-03,
                6.01384044e-05, -2.89286998e-09, 0.00000000e+00]])
