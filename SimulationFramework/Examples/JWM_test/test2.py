import os, sys
import numpy as np
sys.path.append(os.path.abspath(__file__+'/../../../../'))
sys.path.append(os.path.abspath(__file__+'/../../../'))
import SimulationFramework.Framework as Fw

dirname = 'jwm_test2'

Framework = Fw.Framework(dirname, clean = True)
Framework.loadSettings('Lattices/clara400_v12_v3.def')

oldl02phase = Framework.getElement('CLA-L02-CAV', 'phase')
#for n in range(1,5):
#n = 1
for newl02phase in [-10,-5,0,5,10]:
    #Framework.setSubDirectory(dirname + '/l02phase_' + str(newl02phase))
    
    Framework.modifyElement('CLA-L02-CAV', 'phase', oldl02phase + newl02phase)
    Framework.track(startfile='generator', endfile='S07')
