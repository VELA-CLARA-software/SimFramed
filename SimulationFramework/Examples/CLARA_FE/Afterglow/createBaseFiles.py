import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../../'))
from SimulationFramework.Framework import *

''' Run the injector part once if only optimising post-injector parameters'''
def create_base_files(scaling, changes=None):
    framework = Framework('basefiles_'+str(scaling), overwrite=True)
    framework.loadSettings('afterglow.def')
    framework.generator.number_of_particles = 2**(3*scaling)
    if not changes is None:
        framework.load_changes_file(changes)
    framework['S02'].prefix = '../../CLARA/basefiles_' + str(int(scaling)) + '/'
    framework.track(startfile='S02')

## Modify as appropriate! ##
for i in [5]:
    create_base_files(i)
exit()
