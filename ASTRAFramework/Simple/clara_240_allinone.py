from ASTRAInjector import *
import yaml

''' initialise ASTRAFramework and set the base directory to 'example' - Files will be in <working dir>/example/
    If overwrite=False, nothing will be overwritten (useful for testing!) '''
astra = ASTRAInjector('allinone', overwrite=True)

''' Here we *can* set the base value of the ASTRA command
        - defaults to ['astra']
        - in linux we run multi-core using mpiexec
    This *must* be a list    '''
if not os.name == 'nt':
    astra.defineASTRACommand(['mpiexec','-np','6','/opt/ASTRA/astra_MPICH2.sh'])

''' Load a settings file '''
astra.loadSettings('short_240_allinone.settings')
# astra.modifySetting('linac4_field',6156)

''' Create an initial distribution
    - charge is in pC!
    '''
astra.createInitialDistribution(npart=100, charge=250)

''' Apply the new settings '''
astra.applySettings()

''' Run ASTRA '''
astra.runASTRAFiles()
''' The following only runs certain files - this may break if the correct input files do not exist! '''
# astra.runASTRAFiles(files=['test.5'])

''' Run this to create a summary file with all required input files, and the consequent output files'''
# astra.createHDF5Summary(screens=[4929])

''' These commands convert the final bunch to SDDS...'''
ft = feltools('allinone')
sddsfile = ft.convertToSDDS('allinone.in.4929.128')
'''  and then compress it to an RMS value of <dt> '''
# ft.sddsMatchTwiss(sddsfile, 'compressed.sdds', tstdev=5e-13)
