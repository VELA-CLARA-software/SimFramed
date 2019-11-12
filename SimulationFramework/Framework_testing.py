import sys, os
sys.path.append('..')
from SimulationFramework.Framework import Framework

if __name__ == "__main__":
    fw = Framework(directory='test')
    fw.loadSettings('../MasterLattice/Lattices/clara400_v12_FEBE.def')
