# from SAMMcore.Components import *
from SAMMcore.Components import Drift as d
# from SAMMcore.Components import Dipole
from SAMMcore.Components import Quadrupole as Q
from SAMMcore.Components import Screen as S
from SAMMcore.Components import OrbitCorrector as C
from SAMMcore.Components import BeamPositionMonitor as BPM
from SAMMcore.SAMMlab import Beamline
from SAMMcore.SAMMlab import Beam
from SAMMcore.Particles import Electron
from SAMMcore.SAMMlab import PhysicalUnits
import numpy as np

[t1, t2] = np.array([[1], [2]])

beam1 = Beam.Beam(species=Positron.Positron)

beam1.energy = 1.20 * PhysicalUnits.GeV

# help(comp)


# Make VELA libroutine
EBT_S01_DRIFT_01 = d.Drift(length=0.4)
EBT_INJ_MAG_HVCOR_01 = C.OrbitCorrector(field=[0, 0], length=0.05)
EBT_S01_DRIFT_02 = d.Drift(length=0.25)
EBT_INJ_DIA_WCM_01 = d.Drift(length=0.0)
EBT_S01_DRIFT_03 = d.Drift(length=0.15)
EBT_INJ_DIA_BPM_01 = BPM.BeamPositionMonitor(length=0.05)
EBT_S01_DRIFT_04 = d.Drift(length=0.05)
EBT_INJ_MAG_HVCOR_02 = C.OrbitCorrector(field=[0, 0], length=0.05)
EBT_S01_DRIFT_05 = d.Drift(length=0.05)
EBT_INJ_DIA_YAG_01 = S.Screen()
EBT_S01_DRIFT_06 = d.Drift(length=0.185)
EBT_INJ_MAG_QUAD_01 = Q.Quadrupole(length=0.1, gradient=0.3)
EBT_S01_DRIFT_07 = d.Drift(length=0.11)
EBT_INJ_MAG_QUAD_02 = Q.Quadrupole(length=0.1, gradient=-0.1)
EBT_S01_DRIFT_08 = d.Drift(length=0.11)
EBT_INJ_MAG_QUAD_03 = Q.Quadrupole(length=0.1, gradient=0.1)
EBT_S01_DRIFT_09 = d.Drift(length=0.11)
EBT_INJ_MAG_QUAD_04 = Q.Quadrupole(length=0.1, gradient=-0.3)
EBT_S01_DRIFT_10 = d.Drift(length=0.18)
EBT_INJ_DIA_YAG_02 = S.Screen()
EBT_S01_DRIFT_11 = d.Drift(length=0.275)
EBT_INJ_MAG_HVCOR_03 = C.OrbitCorrector(field=[0, 0], length=0.05)
EBT_S01_DRIFT_12 = d.Drift(length=0.05)
EBT_INJ_TDC_01
EBT_S01_DRIFT_13
EBT_INJ_MAG_HVCOR_04






drift1 = Drift.Drift(length=1.5)
quad1 = Quadrupole.Quadrupole(length=0.3, gradient=-0.8)
# dip1 = Dipole.Dipole(field=1, length=0.5, gradient=0.8)


# dip1.theta = 0.2
# dip1.curvature = 3
# print('dip1.curvature = ', dip1.curvature)

# print('dip1.curvature = ', dip1.curvature)

# print dip1.curvature

V1 = Beamline.Beamline(componentlist=[EBT_S01_DRIFT_01,
                                      EBT_INJ_MAG_HVCOR_01,
                                      EBT_S01_DRIFT_02,
                                      EBT_INJ_DIA_WCM_01,
                                      EBT_S01_DRIFT_03,
                                      EBT_INJ_DIA_BPM_01,
                                      EBT_INJ_MAG_HVCOR_02,
                                      EBT_S01_DRIFT_05,
                                      EBT_INJ_DIA_YAG_01,
                                      EBT_S01_DRIFT_06,
                                      EBT_INJ_MAG_QUAD_01,
                                      EBT_S01_DRIFT_07,
                                      EBT_INJ_MAG_QUAD_02,
                                      EBT_S01_DRIFT_08,
                                      EBT_INJ_MAG_QUAD_03,
                                      EBT_S01_DRIFT_09,
                                      EBT_INJ_MAG_QUAD_04,
                                      EBT_S01_DRIFT_10,
                                      EBT_INJ_DIA_YAG_02,
                                      EBT_S01_DRIFT_11,
                                      EBT_INJ_MAG_HVCOR_03,
                                      EBT_S01_DRIFT_12])

# print beamline1.precision
# print beamline1.ComputePositions()

beam1 = Beam.Beam(species=Electron.Electron, energy=4.5 * PhysicalUnits.MeV)
# beam1 = Beam.Beam(species=Positron, energy = 2)

beam1.distance

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
print("beam2 particle1 = ", beam2.particles[0])
print("beam2 particle2 = ", beam2.particles[1])
# print np.array([[9.01351346e-04, -5.97789022e-05, 0.00000000e+00,
#                 0.00000000e+00, -2.50103285e-09, 0.00000000e+00],
#                 [0.00000000e+00, 0.00000000e+00, 1.09921487e-03,
#                 6.01384044e-05, -2.89286998e_09, 0.00000000e+00]])
#
