from ASTRAInjector import *

astra = ASTRAInjector('CSRTest_Base', overwrite=False)
astra.loadSettings('short_240_12b3.settings')
astra.createHDF5Summary(reference='CSRTest_Base')
