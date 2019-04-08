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
index1 = names.index('CLA-S07F-MAG-QUAD-03-K')
index2 = names.index('CLA-S07F-MAG-QUAD-10-K')+1
index3 = names.index('CLA-FEB-MAG-QUAD-13')
# parameter_names = [q for q in names[index1:index2]]
parameter_names.append('CLA-FEB-MAG-QUAD-01')
parameter_names.append('CLA-FEB-MAG-QUAD-02')
parameter_names += [q for q in names[index3:]]
best = parameters
# for n, p in zip(names, parameters)[index3:]:
#     print n, p
# exit()
best = [-0.0744722, 0.47353, -0.0940485, -0.0244515, 0.0257142, -0.10087, 0.57013, -0.61296, #S07
1.55838, -1.09374] #ARC F and D quads
#1.08184, -2.38876, 1.66352, -1.29897, 2.07599, -1.95901, -1.73588, 1.63013, -0.154276, -0.048223, 0.0478511, -0.0479035, 0.0472847]
best = best + parameters[index3:]

with open('transverse_best_changes.yaml', 'r') as infile:
    data = dict(yaml.load(infile, Loader=yaml.UnsafeLoader))
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
        index1 = names.index('CLA-S07F-MAG-QUAD-03-K')
        index2 = names.index('CLA-S07F-MAG-QUAD-10-K') + 1
        index3 = names.index('CLA-FEB-MAG-QUAD-13')
        # self.parameter_names = [[q, 'k1l'] for q in names[index1:index2]]
        self.parameter_names.append(['FODO_F', 'k1l'])
        self.parameter_names.append(['FODO_D', 'k1l'])
        self.parameter_names += [[q, 'k1l'] for q in names[index3:]]
        self.base_files = '../../../CLARA/basefiles_' + str(int(scaling)) + '/'
        self.best_changes = './FEBE_transverse_best_changes.yaml'

    def calculateBeamParameters(self):
        twiss = self.twiss
        self.framework.change_Lattice_Code('All','elegant')
        self.framework['S02'].prefix = self.base_files
        self.framework['S02'].sample_interval = 2**(3*4)
        self.framework.track(startfile='S02')

        constraintsList = {}
        quadkls = self.framework.getElementType('quadrupole','k1l')
        quadlengths = self.framework.getElementType('quadrupole','length')
        quadnames = self.framework.getElementType('quadrupole','objectname')

        twiss.read_elegant_twiss_files( [ self.dirname+'/'+n+'.twi' for n in ['S02', 'L02', 'S03', 'L03', 'S04', 'L4H', 'S05', 'S06', 'L04', 'S07']])
        constraintsListSigmas = {
            'max_xrms': {'type': 'lessthan', 'value': 1e3*twiss['sigma_x'], 'limit': 1, 'weight': 5},
            'max_yrms': {'type': 'lessthan', 'value': 1e3*twiss['sigma_y'], 'limit': 1, 'weight': 5},
            'min_xrms': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_x'], 'limit': 0.1, 'weight': 5},
            'min_yrms': {'type': 'greaterthan', 'value': 1e3*twiss['sigma_y'], 'limit': 0.1, 'weight': 5},
            'last_exn': {'type': 'lessthan', 'value': 1e6*twiss['enx'][-1], 'limit': 0.75, 'weight': 100},
            'last_eyn': {'type': 'lessthan', 'value': 1e6*twiss['eny'][-1], 'limit': 0.75, 'weight': 100},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListSigmas)

        twiss.read_elegant_twiss_files(self.dirname+'/S07.twi')
        tdc_position = self.framework['CLA-S07-TDC-01-R']['position_start'][2]
        tdc_screen_position = self.framework['CLA-S07-DIA-SCR-03-W']['position_start'][2]
        dechirper_position = self.framework['CLA-S07-DCP-01']['position_start'][2]
        constraintsListS07 = {
            # 'tdc_phase_advance': {'type': 'equalto', 'value': twiss.interpolate(tdc_screen_position,'muy') - twiss.interpolate(tdc_position,'muy'), 'limit': 0.25, 'weight': 0.5},
            # 'tdc_screen_beta_y': {'type': 'greaterthan', 'value': twiss.extract_values('beta_y', tdc_position, tdc_screen_position), 'limit': 5, 'weight': 1},
            'dechirper_sigma_x': {'type': 'lessthan', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_x'), 'limit': 0.1, 'weight': 10},
            'dechirper_sigma_y': {'type': 'lessthan', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_y'), 'limit': 0.1, 'weight': 10},
            'dechirper_sigma_xy': {'type': 'equalto', 'value': 1e3*twiss.interpolate(dechirper_position, 'sigma_y') - 1e3*twiss.interpolate(dechirper_position, 'sigma_x'), 'limit': 0.0, 'weight': 10},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListS07)

        twiss.read_elegant_twiss_files( self.dirname+'/S07F.twi' )
        constraintsListFEBEStart = {
            'betax': {'type': 'equalto', 'value': twiss['beta_x'][-1], 'limit': 0.706961, 'weight': 25},
            'betay': {'type': 'equalto', 'value': twiss['beta_y'][-1], 'limit': 3.94176, 'weight': 25},
            'alphax': {'type': 'equalto', 'value': twiss['alpha_x'][-1], 'limit': -1.40098, 'weight': 25},
            'alphay': {'type': 'equalto', 'value': twiss['alpha_y'][-1], 'limit': 0.872449, 'weight': 25},
            'etax': {'type': 'equalto', 'value': twiss.elegant['etax'][-1], 'limit': 0., 'weight': 500},
            'etaxp': {'type': 'equalto', 'value': twiss.elegant['etaxp'][-1], 'limit': 0., 'weight': 500},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListFEBEStart)

        twiss.read_elegant_twiss_files( self.dirname+'/FEBE.twi' )
        ipindex = list(twiss.elegant['ElementName']).index('CLA-FEB-W-FOCUS-01')
        constraintsListFEBE = {
            'max_betax': {'type': 'lessthan', 'value': twiss['beta_x'], 'limit': 30, 'weight': 5},
            'max_betay': {'type': 'lessthan', 'value': twiss['beta_y'], 'limit': 30, 'weight': 5},
            'ip_sigmax': {'type': 'lessthan', 'value': 1e3*twiss.elegant['Sx'][ipindex], 'limit': 0.05, 'weight': 5},
            'ip_sigmay': {'type': 'lessthan', 'value': 1e3*twiss.elegant['Sy'][ipindex], 'limit': 0.05, 'weight': 5},
            'ip_alphax': {'type': 'equalto', 'value': twiss.elegant['alphax'][ipindex], 'limit': 0., 'weight': 5},
            'ip_alphay': {'type': 'equalto', 'value': twiss.elegant['alphay'][ipindex], 'limit': 0., 'weight': 5},
            'ip_etax': {'type': 'equalto', 'value': twiss.elegant['etax'][ipindex], 'limit': 0., 'weight': 5000},
            'ip_etaxp': {'type': 'equalto', 'value': twiss.elegant['etaxp'][ipindex], 'limit': 0., 'weight': 5000},
            'ip_enx': {'type': 'lessthan', 'value': 1e6*twiss.elegant['enx'][ipindex], 'limit': 2, 'weight': 15},
            'ip_eny': {'type': 'lessthan', 'value': 1e6*twiss.elegant['eny'][ipindex], 'limit': 0.5, 'weight': 25},
        }
        constraintsList = merge_two_dicts(constraintsList, constraintsListFEBE)

        fitness = self.cons.constraints(constraintsList)
        if self.verbose:
            print self.cons.constraintsList(constraintsList)
        return fitness


if __name__ == "__main__":
        fit = FEBE_Transverse('./FEBE.def', scaling=6)
        fit.setChangesFile(['./nelder_mead_best_changes.yaml', './transverse_best_changes.yaml'])
        fit.verbose = False
        fit.Nelder_Mead(best, step=0.1)
        # fit.Simplex(best)
