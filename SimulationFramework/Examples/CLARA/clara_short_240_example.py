import os
import sys
sys.path.append(os.path.abspath(__file__+'/../../../../'))
from SimulationFramework.Framework import *

ncpu = 20

framework = Framework('Short-240_Combined', clean=True)
framework.loadSettings('./Short-240_Combined.def')
framework.load_changes_file(['./Elegant_Genesis/Short-240/nelder_mead_best_changes.yaml', './Elegant_Genesis/Short-240/transverse_best_changes.yaml'])
framework.change_Lattice_Code('All','elegant')
# framework['S02'].prefix = '../basefiles_6/'
framework['S02'].sample_interval = 2**(1*3)
framework.save_changes_file(filename=framework.subdirectory+'/changes.yaml')
framework.track(startfile='S02', endfile='S07')

framework = Framework('Short-240_Separated', clean=True)
framework.loadSettings('./Short-240_Separated.def')
framework.load_changes_file(['./Elegant_Genesis/Short-240/nelder_mead_best_changes.yaml', './Elegant_Genesis/Short-240/transverse_best_changes.yaml'])
framework.change_Lattice_Code('All','elegant')
# framework['S02'].prefix = '../basefiles_6/'
framework['S02'].sample_interval = 2**(1*3)
framework.save_changes_file(filename=framework.subdirectory+'/changes.yaml')
framework.track(startfile='S02', endfile='S07')
