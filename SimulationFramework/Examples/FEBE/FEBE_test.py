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
        self.sample_interval = 2**(3*4)
        self.base_files = '../../CLARA/basefiles_'+str(self.scaling)+'/'
        self.clean = True
        self.doTracking = True
        self.change_to_elegant = True
        self.change_to_astra = False
        self.change_to_gpt = False
        # self.start_lattice = 'FEBE'

    def before_tracking(self):
            elements = self.framework.elementObjects.values()
            for e in elements:
                e.lsc_enable = True
                e.lsc_bins = 100
                e.current_bins = 0
                e.longitudinal_wakefield_enable = True
                e.transverse_wakefield_enable = True
                e.smoothing_half_width = 2
                pass
            lattices = self.framework.latticeObjects.values()
            for l in lattices:
                l.lscDrifts = True
                l.lsc_bins = 100
                l.lsc_high_frequency_cutoff_start = 0.25
                l.lsc_high_frequency_cutoff_end = 0.33
                pass
            self.framework['FEBE'].betax = 0.74306
            self.framework['FEBE'].betay = 3.96111
            self.framework['FEBE'].alphax = -0.623844
            self.framework['FEBE'].alphay = 0.872959

            BPM1_replacement_index = list(self.framework['FEBE'].elements.keys()).index('CLA-FEB-DIA-BPM-01')
            self.framework.add_Element(name='CLA-FEBE-TDC-01', type='rf_deflecting_cavity', **{'length': self.framework['CLA-FEB-DIA-BPM-01'].length, 'phase': 0.0, 'global_rotation': [0, 0, 0],
            'position_end': self.framework['CLA-FEB-DIA-BPM-01'].position_end, 'position_start': self.framework['CLA-FEB-DIA-BPM-01'].position_start})
            self.framework['FEBE'].insert_element()
            print(self.framework['CLA-FEBE-TDC-01'])
            exit()
            # self.framework['laser-heater'].angle = 0
            # self.framework['CLA-L4H-CAV'].field_amplitude = 0
            # self.framework['bunch_compressor'].angle = 1.18*0.12850503205442976

if __name__ == "__main__":
    opt = FEBE()
    opt.set_changes_file(['./transverse_best_changes_upto_S07.yaml', './S07_transverse_best_changes.yaml'])
    opt.set_lattice_file('./FEBE_Single.def')
    opt.set_start_file('PreFEBE')
    opt.load_best('./nelder_mead_best_changes.yaml')
    opt.Example(dir='example_GPT')
    # opt.Simplex()
