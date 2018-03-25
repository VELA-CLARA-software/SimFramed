import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname( os.path.abspath(__file__)))))
from FrameworkTest.framework import *

lattice = framework('CLARA')
lattice.loadSettings('Lattices/clara400_v12_elegant.def')
# print lattice.elements
print lattice.getElement('CLA-S02-MAG-QUAD-01').write_ASTRA(1)
for i, q in enumerate(lattice.quadrupoles,1):
    print q.write_ASTRA(i)
