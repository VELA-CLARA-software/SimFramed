from ASTRAGeneral import *
import numpy as np

npart=2**(3*4)

grids = getGrids(npart=npart)
print grids.getGridSizes()


astra = ASTRA('1')
if not os.name == 'nt':
    astra.defineASTRACommand(['mpiexec','-np','20','/opt/ASTRA/astra_MPICH2.sh'])
astra.loadSettings('clara_v11.def')
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
astra.createInitialDistribution(npart=npart,charge=250)
astra.createASTRAFiles()
astra.runASTRAFiles()
# astra.runASTRAFiles(['VBC','S06','L04','S07'])
