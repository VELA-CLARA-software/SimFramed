import numpy as np
import os, sys
sys.path.append(os.path.abspath(__file__+'/../../../../../'))
from FitnessFunc_Transverse import *
import operator
import random
import csv
from scipy.optimize import minimize
import shutil
import uuid
from functools import partial


def optfunc(args, dir=None, **kwargs):
    if dir == None:
        with TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            fit = fitnessFunc(args, tmpdir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    else:
            fit = fitnessFunc(args, dir, **kwargs)
            fitvalue = fit.calculateBeamParameters()
    return (fitvalue,)


framework = Framework('transverse_best', overwrite=False)
framework.loadSettings('Lattices/clara400_v12_v3.def')
best = framework.getElementType('quadrupole','k1')
#
# allbest=[]
# with open('transverse_best_Short_240_solutions.csv.tmp','r') as infile:
#     reader = csv.reader(infile, quoting=csv.QUOTE_NONE, skipinitialspace=True)
#     for row in reader:
#         allbest.append(row)
# best = map(lambda x: float(x), allbest[0])

print 'starting values = ', best

# print optfunc(best, dir=os.getcwd()+'/test_transverse', scaling=6, overwrite=True, verbose=True, summary=False)
# exit()


optfunc3 = partial(optfunc, dir=None, scaling=3, post_injector=True, parameterfile='clara_longitudinal_best.yaml', verbose=True)
optfunc4 = partial(optfunc, dir=None, scaling=4, post_injector=True, parameterfile='clara_longitudinal_best.yaml')

print optfunc3(best)
exit()

with open('best_solutions_running_simplex.csv','w') as out:
    pass
res = minimize(optfunc3, best, method='nelder-mead', options={'xtol': 1e-8, 'disp': True, 'adaptive': True})
print res

try:
    print 'best fitness = ', optfunc(res.x, dir=os.getcwd()+'/CLARA_best_longitudinal_simplex', scaling=6, overwrite=True, verbose=True, summary=True, post_injector=False)
except:
    pass

saveParameterFile(res.x)
