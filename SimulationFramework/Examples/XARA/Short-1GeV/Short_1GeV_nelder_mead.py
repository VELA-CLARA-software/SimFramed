import sys, os
import numpy as np
sys.path.append('./../')
from Optimise_Genesis_Elegant import Optimise_Genesis_Elegant
import ruamel.yaml as yaml

class Short1GeV(Optimise_Genesis_Elegant):

    parameter_names = [
        # ['CLA-L01-CAV', 'field_amplitude'],
        # ['CLA-L01-CAV', 'phase'],
        ['CLA-L02-CAV', 'field_amplitude'],
        ['CLA-L02-CAV', 'phase'],
        ['CLA-L03-CAV', 'field_amplitude'],
        ['CLA-L03-CAV', 'phase'],
        ['CLA-L4H-CAV', 'field_amplitude'],
        ['CLA-L4H-CAV', 'phase'],
        ['CLA-L04-CAV-01', 'field_amplitude'],
        ['CLA-L04-CAV-01', 'phase'],
        ['CLA-L04-CAV-02', 'field_amplitude'],
        ['CLA-L04-CAV-02', 'phase'],
        ['CLA-L04-CAV-03', 'field_amplitude'],
        ['CLA-L04-CAV-03', 'phase'],
        ['CLA-L04-CAV-04', 'field_amplitude'],
        ['CLA-L04-CAV-04', 'phase'],
        ['bunch_compressor', 'angle'],
        ['CLA-S07-DCP-01', 'factor'],
    ]

    def __init__(self):
        super(Short1GeV, self).__init__()
        self.CLARA_dir = os.path.relpath(__file__+'/../../../../CLARA/')
        self.scaling = 6
        self.sample_interval=2**(3*1)
        self.base_files = '../../../../CLARA/basefiles_6/'
        self.genesis_file = 'xara_4nm_td.in'
        self.verbose = False
        self.alphax = -0.08050797
        self.betax = 4.58986
        self.alphay = -0.05371444
        self.betay = 2.52698

    def load_best(self, filename):
        with open(filename, 'r') as infile:
            data = dict(yaml.load(infile, Loader=yaml.UnsafeLoader))
            # best = [data[n]['k1l'] for n in parameter_names]
            best = []
            for n, p in self.parameter_names:
                if n in data:
                    best.append(data[n][p])
                elif n == 'bunch_compressor' and p == 'set_angle':
                    best.append(data['CLA-VBC-MAG-DIP-01']['angle'])
                else:
                    print(n, p)
                    if not hasattr(self, 'framework'):
                        self.framework = fw.Framework(None)
                        self.framework.loadSettings(self.lattice)
                    best.append(self.framework[n][p])
            self.best = best
        return best

    def calculate_constraints(self):

        if len(self.linac_fields) > 0:
            slinac_fields = 1e-6*np.array([self.linac_fields[i] for i in [0,1]])
            x4hlinac_fields = 1e-6*np.array([self.linac_fields[i] for i in [2]])
            xlinac_fields = 1e-6*np.array([self.linac_fields[i] for i in [4,5,6]])
        else:
            slinac_fields = []
            x4hlinac_fields = []
            xlinac_fields = []

        r = self.resultsDict
        pos12m = list(r['g'].z).index(10.08)
        pos20m = list(r['g'].z).index(16.08)
        print 'Momentum = ', r['momentum']
        constraintsList = {
            'brightness': {'type': 'greaterthan', 'value': r['brightness'], 'limit': 0.15, 'weight': 0},
            'bandwidth': {'type': 'lessthan', 'value': r['b'], 'limit': 0.1, 'weight': 3},
            'pulse_energy': {'type': 'greaterthan', 'value': 1e2*r['e'], 'limit': 150, 'weight': 4},
            'bandwidth_end': {'type': 'lessthan', 'value': 1e2*abs(r['g'].spectrum_lamdwidth_std[pos20m]), 'limit': 0.25, 'weight': 1},
            'pulse_energy_end': {'type': 'greaterthan', 'value': 1e6*abs(r['g'].energy[pos20m]), 'limit': 250, 'weight': 2},
            'max_brightness_position': {'type': 'lessthan', 'value': abs(r['l']), 'limit': 10, 'weight': 2.5},
            'min_brightness_position': {'type': 'greaterthan', 'value': abs(r['l']), 'limit': 7, 'weight': 0.5},
            'energy_min': {'type': 'greaterthan', 'value': abs(r['momentum']), 'limit': 950, 'weight': 5},
            'energy_max': {'type': 'lessthan', 'value': abs(r['momentum']), 'limit': 1050, 'weight': 5},
            'field_max_s': {'type': 'lessthan', 'value': slinac_fields, 'limit': 32, 'weight': 3},
            'field_max_x': {'type': 'lessthan', 'value': xlinac_fields, 'limit': 80, 'weight': 3},
            'field_max_x4h': {'type': 'lessthan', 'value': x4hlinac_fields, 'limit': 35, 'weight': 3},
        }
        return constraintsList

if __name__ == "__main__":
    opt = Short1GeV()
    opt.set_changes_file(['../nelder_mead_best_changes.yaml', '../transverse_best_changes.yaml'])
    opt.set_lattice_file('Lattices/claraX400_v12_80MVm_Elegant.def')
    opt.load_best('../nelder_mead_best_changes.yaml')
    opt.start_lattice = 'CLARAX'
    opt.runElegant = True
    opt.runGenesis = True
    opt.Nelder_Mead(step=[5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5, 5e6, 5,  5e6, 5,  5e6, 5, 0.005, 0.1])
