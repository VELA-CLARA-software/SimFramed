import sys, os
sys.path.append('../../../')
from SimulationFramework.ClassFiles.Optimise_longitudinal_Elegant import Optimise_Elegant
import SimulationFramework.Modules.read_twiss_file as rtf
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Example Tracking')
parser.add_argument('type', default='nelder_mead')

class FEBE(Optimise_Elegant):

    # injector_startingvalues = [-9.,0.345,2.1e7,-16.,0.052500000000000005,-0.05]
    # startingvalues = best = np.array([ 3.13845650e+07, -2.33062481e+01,  2.96752546e+07, -3.41502595e+00,
    #     2.74842883e+07,  1.87482967e+02,  3.11918859e+07,  5.07160187e+01,
    #    -1.22393267e-01,  6.12784140e-01])

    parameter_names = [
        # ['startcharge','charge'],
        ['CLA-L01-CAV', 'field_amplitude'],
        ['CLA-L01-CAV', 'phase'],
        ['CLA-L02-CAV', 'field_amplitude'],
        ['CLA-L02-CAV', 'phase'],
        ['CLA-L03-CAV', 'field_amplitude'],
        ['CLA-L03-CAV', 'phase'],
        ['CLA-L4H-CAV', 'field_amplitude'],
        ['CLA-L4H-CAV', 'phase'],
        ['CLA-L04-CAV', 'field_amplitude'],
        ['CLA-L04-CAV', 'phase'],
        ['bunch_compressor', 'angle'],
        ['CLA-S07-DCP-01', 'factor'],
        # ['FODO_D', 'k1l'],
        # ['FODO_F', 'k1l'],
    ]

    def __init__(self):
        super(FEBE, self).__init__()
        self.parameter_names.append(['FODO_F', 'k1l'])
        self.parameter_names.append(['FODO_D', 'k1l'])
        self.scaling = 6
        self.sample_interval = 2**(3*0)
        self.base_files = '../basefiles_'+str(self.scaling)+'/'
        self.clean = True
        self.doTracking = True
        self.change_to_elegant = True
        self.change_to_astra = False
        self.change_to_gpt = False

    def before_tracking(self):
            self.framework.defineElegantCommand(scaling=4)
            csrbins = int(round(2**(3*self.scaling) / self.sample_interval / 128, -1))
            lscbins = int(round(2**(3*self.scaling) / self.sample_interval / 256, -1))
            elements = self.framework.elementObjects.values()
            for e in elements:
                e.lsc_enable = True
                e.lsc_bins = lscbins
                e.current_bins = 0
                e.csr_bins = csrbins
                e.longitudinal_wakefield_enable = True
                e.transverse_wakefield_enable = True
                e.smoothing_half_width = 1
                e.lsc_high_frequency_cutoff_start = 0.2
                e.lsc_high_frequency_cutoff_end = 0.25
                e.smoothing = 1
                # e.end1_focus = 0
                # e.end2_focus = 0
                # e.body_focus_model = "None"
            lattices = self.framework.latticeObjects.values()
            for l in lattices:
                l.lscDrifts = True
                l.lsc_bins = lscbins
                l.lsc_high_frequency_cutoff_start = 0.2
                l.lsc_high_frequency_cutoff_end = 0.25
                l.smoothing_half_width = 1
                l.smoothing = 1
            # self.framework['FEBE'].betax = 0.748377
            # self.framework['FEBE'].betay = 3.96107
            # self.framework['FEBE'].alphax = -0.61447
            # self.framework['FEBE'].alphay = 0.872954

            # BPM1_replacement_index = list(self.framework['FEBE'].elements.keys()).index('CLA-FEB-DIA-BPM-01')
            # self.framework.add_Element(name='CLA-FEBE-TDC-01', type='rf_deflecting_cavity', **{'length': self.framework['CLA-FEB-DIA-BPM-01'].length, 'phase': 0.0, 'global_rotation': [0, 0, 0],
            # 'position_end': self.framework['CLA-FEB-DIA-BPM-01'].position_end, 'position_start': self.framework['CLA-FEB-DIA-BPM-01'].position_start})
            # self.framework['FEBE'].insert_element()
            # print(self.framework['CLA-FEBE-TDC-01'])
            # exit()
            # self.framework['laser-heater'].angle = 0
            # self.framework['CLA-L4H-CAV'].field_amplitude = 0
            # self.framework['bunch_compressor'].angle = 1.18*0.12850503205442976

if __name__ == "__main__":
    opt = FEBE()
    opt.set_changes_file(['./transverse_best_changes_upto_S07_250pC.yaml', './S07_transverse_best_changes_250pC.yaml','./FEBE_transverse_best_changes.yaml'])
    opt.set_lattice_file('./FEBE_Single_L01.def')
    opt.set_start_file('PreFEBE')
    args = parser.parse_args()
    if args.type == 'nelder_mead':
        print('Using a Nelder-Mead Track...')
        opt.load_best('./nelder_mead_best_changes_250pC.yaml')
    elif args.type == 'simplex':
        print('Using a SciPy Simplex Track...')
        opt.load_best('./simplex_best_changes_250pC.yaml')
    opt.Example(dir='example_250pC')
    opt.framework.save_lattice(directory='example_250pC/')
