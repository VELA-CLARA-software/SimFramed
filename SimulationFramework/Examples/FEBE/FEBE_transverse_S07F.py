import os, sys
sys.path.append('../../../')
import SimulationFramework.Framework as fw
from SimulationFramework.Modules.nelder_mead import nelder_mead
from SimulationFramework.Examples.CLARA.Elegant.Optimise_transverse import Optimise_transverse
from SimulationFramework.Modules.merge_two_dicts import merge_two_dicts
from ruamel import yaml

framework = fw.Framework(None)
framework.loadSettings('FEBE.def')
parameters = framework.getElementType('quadrupole','k1l')
names = framework.getElementType('quadrupole','objectname')
index1 = names.index('CLA-S07-MAG-QUAD-01')
index2 = names.index('CLA-S07-MAG-QUAD-10')+1
index3 = names.index('CLA-FEB-MAG-QUAD-13')
parameter_names = [q for q in names[index1:index2]]
# parameter_names.append('CLA-FEB-MAG-QUAD-01')
# parameter_names.append('CLA-FEB-MAG-QUAD-02')
# parameter_names += [q for q in names[index3:]]
best = parameters
print(parameter_names)
# for n, p in zip(names, parameters)[index3:]:
#     print n, p
# exit()
# best = [-0.0744722, 0.47353, -0.0940485, -0.0244515, 0.0257142, -0.10087, 0.57013, -0.61296, #S07
# 1.55838, -1.09374] #ARC F and D quads
# #1.08184, -2.38876, 1.66352, -1.29897, 2.07599, -1.95901, -1.73588, 1.63013, -0.154276, -0.048223, 0.0478511, -0.0479035, 0.0472847]
# best = best + parameters[index3:]
#
with open('S07_transverse_best_changes.yaml', 'r') as infile:
    data = dict(yaml.load(infile, Loader=yaml.UnsafeLoader))
    # data=[]
    # best = [data[n]['k1l'] for n in parameter_names]
    best = []
    for n in parameter_names:
        if n in data:
            best.append(data[n]['k1l'])
        else:
            best.append(framework[n]['k1l'])

class FEBE_Transverse(Optimise_transverse):

    def __init__(self, lattice='FEBE.def', scaling=6):
        super(FEBE_Transverse, self).__init__(lattice=lattice, scaling=scaling)
        names = framework.getElementType('quadrupole','objectname')
        index1 = names.index('CLA-S07-MAG-QUAD-01')
        index2 = names.index('CLA-S07-MAG-QUAD-10') + 1
        index3 = names.index('CLA-FEB-MAG-QUAD-13')
        self.parameter_names = []
        self.parameters = []
        self.parameter_names = [q for q in names[index1:index2]]
        self.parameters = [[q, 'k1l'] for q in names[index1:index2]]
        # self.parameter_names.append('CLA-FEB-MAG-QUAD-01')
        # self.parameters.append(['FODO_F', 'k1l'])
        # self.parameter_names.append('CLA-FEB-MAG-QUAD-02')
        # self.parameters.append(['FODO_D', 'k1l'])
        # self.parameter_names += [q for q in names[index3:]]
        # self.parameters += [[q, 'k1l'] for q in names[index3:]]
        self.base_files = '../../../CLARA/basefiles_' + str(int(scaling)) + '/'
        self.best_changes = './S07_transverse_best_changes.yaml'
        self.start_file = 'PreFEBE'

    def calculateBeamParameters(self):
        twiss = self.twiss
        self.framework.change_Lattice_Code('All','elegant')
        self.framework[self.start_file].prefix = self.base_files
        self.framework[self.start_file].sample_interval = 2**(3*4)
        self.framework.track(startfile=self.start_file)

        constraintsList = {}
        quadkls = self.framework.getElementType('quadrupole','k1l')
        quadlengths = self.framework.getElementType('quadrupole','length')
        quadnames = self.framework.getElementType('quadrupole','objectname')

        # twiss.read_elegant_twiss_files( [ self.dirname+'/'+n+'.twi' for n in ['S07', 'S07F']])
        twiss.read_elegant_twiss_files( self.dirname+'/'+self.start_file+'.twi' )
        constraintsListSigmas = {
            'max_xrms': {'type': 'lessthan', 'value': max(1e3*twiss['sigma_x']), 'limit': 1, 'weight': 5},
            'max_yrms': {'type': 'lessthan', 'value': max(1e3*twiss['sigma_y']), 'limit': 1, 'weight': 5},
            'min_xrms': {'type': 'greaterthan', 'value': min(1e3*twiss['sigma_x']), 'limit': 0.05, 'weight': 50},
            'min_yrms': {'type': 'greaterthan', 'value': min(1e3*twiss['sigma_y']), 'limit': 0.05, 'weight': 50},
            'last_exn': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 2, 'weight': 0},
            'last_eyn': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.75, 'weight': 0},
            'max_betax': {'type': 'lessthan', 'value': max(twiss['beta_x']), 'limit': 250, 'weight': 5},
            'max_betay': {'type': 'lessthan', 'value': max(twiss['beta_y']), 'limit': 250, 'weight': 5},
            # 'max_alphax': {'type': 'lessthan', 'value': max(abs(twiss['alpha_x'])), 'limit': 10, 'weight': 5},
            # 'max_alphay': {'type': 'lessthan', 'value': max(abs(twiss['alpha_y'])), 'limit': 10, 'weight': 5},
        }
        # print self.cons.constraintsList(constraintsListSigmas)
        constraintsList = merge_two_dicts(constraintsList, constraintsListSigmas)

        # twiss.read_elegant_twiss_files(self.dirname+'/S07F.twi')
        dechirper_position = self.framework['CLA-S07-DCP-01']['position_start'][2]
        dcpindex = list(twiss.elegant['ElementName']).index('CLA-S07-DCP-01')
        constraintsListS07 = {
            'dechirper_sigma_x': {'type': 'lessthan', 'value': 1e3*twiss.elegant['Sx'][dcpindex], 'limit': 0.1, 'weight': 2},
            'dechirper_sigma_y': {'type': 'lessthan', 'value': 1e3*twiss.elegant['Sy'][dcpindex], 'limit': 0.1, 'weight': 2},
            'dechirper_sigma_xy': {'type': 'equalto', 'value': 1e3*twiss.elegant['Sy'][dcpindex] - 1e3*twiss.elegant['Sx'][dcpindex], 'limit': 0.0, 'weight': 2},
        }
        # constraintsList = merge_two_dicts(constraintsList, constraintsListS07)

        # twiss.read_elegant_twiss_files( self.dirname+'/S07F.twi' )
        constraintsListFEBEStart = {
            'betax': {'type': 'equalto', 'value': twiss['beta_x'][-1], 'limit': 0.706961, 'weight': 50},
            'betay': {'type': 'equalto', 'value': twiss['beta_y'][-1], 'limit': 3.94176, 'weight': 50},
            'alphax': {'type': 'equalto', 'value': twiss['alpha_x'][-1], 'limit': -1.40098, 'weight': 50},
            'alphay': {'type': 'equalto', 'value': twiss['alpha_y'][-1], 'limit': 0.872449, 'weight': 50},
            # 'etax': {'type': 'equalto', 'value': twiss.elegant['etax'][-1], 'limit': 0., 'weight': 500},
            # 'etaxp': {'type': 'equalto', 'value': twiss.elegant['etaxp'][-1], 'limit': 0., 'weight': 500},
        }
        # print self.cons.constraintsList(constraintsListFEBEStart)
        constraintsList = merge_two_dicts(constraintsList, constraintsListFEBEStart)

        self.constraintsList = constraintsList
        fitness = self.cons.constraints(constraintsList)
        if self.verbose:
            print(self.cons.constraintsList(constraintsList))
        return fitness


if __name__ == "__main__":
        fit = FEBE_Transverse('./FEBE_Single.def', scaling=6)
        fit.setChangesFile(['./nelder_mead_best_changes.yaml','./transverse_best_changes_upto_S07.yaml'])
        fit.verbose = False
        fit.Nelder_Mead(best, step=0.01)
        # fit.Simplex(best)
