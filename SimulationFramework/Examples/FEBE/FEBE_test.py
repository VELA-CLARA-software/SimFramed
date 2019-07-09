import sys, os
sys.path.append('../../../')
from SimulationFramework.ClassFiles.Optimise_longitudinal_Elegant import Optimise_Elegant
import SimulationFramework.Modules.read_twiss_file as rtf
import numpy as np

class FEBE(Optimise_Elegant):

    # injector_startingvalues = [-9.,0.345,2.1e7,-16.,0.052500000000000005,-0.05]
    # startingvalues = best = np.array([ 3.13845650e+07, -2.33062481e+01,  2.96752546e+07, -3.41502595e+00,
    #     2.74842883e+07,  1.87482967e+02,  3.11918859e+07,  5.07160187e+01,
    #    -1.22393267e-01,  6.12784140e-01])

    def __init__(self):
        super(FEBE, self).__init__()
        self.parameter_names.append(['FODO_F', 'k1l'])
        self.parameter_names.append(['FODO_D', 'k1l'])
        self.scaling = 6
        self.sample_interval = 2**(3*1)
        self.base_files = '../../CLARA/basefiles_'+str(self.scaling)+'/'
        self.clean = True
        self.doTracking = True

if __name__ == "__main__":
    opt = FEBE()
    opt.set_changes_file(['./transverse_best_changes_upto_S07.yaml', './S07_transverse_best_changes.yaml', './FEBE_transverse_best_changes.yaml'])
    opt.set_lattice_file('./FEBE_Single.def')
    opt.set_start_file('PreFEBE')
    opt.load_best('./nelder_mead/iteration_201/changes.yaml')
    opt.Example()
    # opt.Simplex()
