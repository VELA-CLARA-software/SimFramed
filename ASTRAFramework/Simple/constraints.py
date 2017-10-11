import numpy as np

class constraintsClass():

    def lessthan(self, type, value, limit, weight=1):
        if value > limit:
            return (weight*np.abs(value - limit))**2
        else:
            return 0

    def greaterthan(self, type, value, limit, weight=1):
        if value < limit:
            return (weight*np.abs(value - limit))**2
        else:
            return 0

    def equalto(self, type, value, limit, weight=1):
        return (weight*np.abs(value - limit))**2

    def constraints(self, constraints={}):
        ans = 0
        for k, v in constraints.items():
            if hasattr(self, v['type'].lower()):
                ans += getattr(self, v['type'].lower())(**v)
        return np.sqrt(ans)
