import os, sys
import numpy as np
sys.path.append(os.path.abspath(__file__+'/../../../../'))
import SimulationFramework.Framework as Fw

Framework = Fw.Framework('jwm_test1')
Framework.loadSettings('Lattices/clara400_v12_v3.def')
Framework.track(startfile='generator', endfile='S07')
