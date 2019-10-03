import os, sys
sys.path.append('../../../../')
import SimulationFramework.Framework as fw
from SimulationFramework.Modules.nelder_mead import nelder_mead
from SimulationFramework.ClassFiles.Optimise_transverse import Optimise_transverse
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from ruamel import yaml

framework = fw.Framework(None)
framework.loadSettings('afterglow.def')
parameters = framework.getElementType('quadrupole','k1l')
names = framework.getElementType('quadrupole','objectname')
index1 = names.index('EBT-BA1-MAG-QUAD-04')
parameter_names = []
# parameter_names = [q for q in names[index1:index2]]
# parameter_names.append('FODO_F')
# parameter_names.append('FODO_D')
parameter_names += [q for q in names[index1:]]
best = [1,-1,1,-4]

with open('transverse_best_changes.yaml', 'r') as infile:
    data = dict(yaml.load(infile, Loader=yaml.UnsafeLoader))
    # best = [data[n]['k1l'] for n in parameter_names]
    best = []
    for n in parameter_names:
        if n in data:
            best.append(data[n]['k1l'])
        else:
            print(n)
            best.append(framework[n]['k1l'])

class Afterglow_Transverse(Optimise_transverse):

    def __init__(self, lattice='afterglow.def', scaling=6):
        super(Afterglow_Transverse, self).__init__(lattice=lattice, scaling=scaling)
        names = framework.getElementType('quadrupole','objectname')
        self.parameter_names = []
        self.parameters = []
        self.parameter_names += [q for q in names[index1:]]
        self.parameters += [[q, 'k1l'] for q in names[index1:]]
        self.base_files = '../../'
        self.best_changes = './transverse_best_changes.yaml'
        # self.start_file = 'S02'
        self.start_file = 'PostBA1'

    def calculateBeamParameters(self):
        twiss = self.twiss
        self.framework.change_Lattice_Code('All','elegant')
        self.framework.defineElegantCommand(location=['elegant'])
        self.framework[self.start_file].prefix = self.base_files
        self.framework[self.start_file].sample_interval = 2**(3*4)
        self.framework.track(startfile=self.start_file)

        constraintsList = {}

        twiss.reset_dicts()

        for lat in ['PostBA1']:
            quadkls = self.framework[lat].getElementType('quadrupole','k1l')
            quadlengths = self.framework[lat].getElementType('quadrupole','length')
            constraintsListQuads = {
                'max_k_'+lat: {'type': 'lessthan', 'value': [abs(k) for k, l in zip(quadkls, quadlengths)], 'limit': 5.0, 'weight': 75},

            }
            constraintsList = merge_two_dicts(constraintsList, constraintsListQuads)

        twiss.read_elegant_twiss_files( [ self.dirname+'/PostBA1.twi'])
        constraintsListTwiss = {
            'max_betax': {'type': 'lessthan', 'value': max(twiss['beta_x']), 'limit': 10, 'weight': 150},
            'max_betay': {'type': 'lessthan', 'value': max(twiss['beta_y']), 'limit': 10, 'weight': 150},
            'dump_etax': {'type': 'equalto', 'value': twiss['eta_x'][-1], 'limit': 0.67, 'weight': 5000},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListTwiss)

        self.constraintsList = constraintsList
        fitness = self.cons.constraints(constraintsList)
        if self.verbose:
            print(self.cons.constraintsList(constraintsList))
        return fitness


if __name__ == "__main__":
        fit = Afterglow_Transverse('./afterglow.def', scaling=6)
        # fit.setChangesFile([])
        fit.verbose = False
        fit.Nelder_Mead(best, step=0.1)
        # fit.Simplex(best)
