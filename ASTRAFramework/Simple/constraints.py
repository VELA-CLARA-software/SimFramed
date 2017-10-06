import numpy as np

def lessThan(value, limit, weight=1):
    if value > limit:
        return (weight*np.abs(value - limit))**2
    else:
        return 0

def greaterThan(value, limit, weight=1):
    if value < limit:
        return (weight*np.abs(value - limit))**2
    else:
        return 0

def equalTo(value, limit, weight=1):
    return (weight*np.abs(value - limit))**2

def constraints(constraints={}):
    ans = []
    for k, v in constraints.items():
        ans.append(locals()[k](*v))
    map(lambda x: x**2, ans)
