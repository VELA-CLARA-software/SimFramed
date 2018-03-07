from Framework import *
import numpy as np

npart=2**(3*5)

grids = getGrids(npart=npart)
# print grids.getGridSizes()


framework = Framework('2', overwrite=True)
if not os.name == 'nt':
    framework.astra.defineASTRACommand(['mpiexec','-np','20','/opt/ASTRA/astra_MPICH2.sh'])
framework.loadSettings('clara400_v12.def')
# astra.elements['CLA-S02-MAG-QUAD-01']['K'] = input
# print astra.createASTRAChicane('laser-heater',1)

# print astra.createASTRAQuad('QUAD-01',1)

# print astra.createASTRASolenoid('SOL-01',1)

# print astra.createASTRACavity('GUN10',1)
# print astra.createASTRACavity('LINAC-01',2)

# print astra.createASTRAChargeBlock()

# print astra.createASTRANewRunBlock()
# print astra.createASTRAChargeBlock()
# print astra.createASTRAScanBlock()
# print astra.createASTRAApertureBlock()
# print astra.createASTRACavityBlock()
# print astra.createASTRASolenoidBlock()
# print astra.createASTRAQuadrupoleBlock()
framework.astra.createInitialDistribution(npart=npart,charge=250)
gun = framework.getElement('CLA-HRG1-GUN-CAV')
gun['phase'] = -5
framework.createInputFiles()
framework.runInputFiles(['injector400'])
# astra.runASTRAFiles(['VBC','S06','L04','S07'])
