from constraints import *

cons = constraintsClass()

constraintsList = {'lessThan':{'value': 23, 'limit': 0, 'weight': 1},
'greaterThan': {'value':0,'limit':12,'weight':1}}

print cons.constraints(constraintsList)
