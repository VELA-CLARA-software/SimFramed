from collections import OrderedDict
from copy import copy

def merge_two_dicts(y, x):
    '''Combine to dictionaries: first dictionary overwrites keys in the second dictionary'''
    if not isinstance(x, (dict, OrderedDict)) and not isinstance(y, (dict, OrderedDict)):
        return OrderedDict()
    elif not isinstance(x, (dict, OrderedDict)):
        return y
    elif not isinstance(y, (dict, OrderedDict)):
        return x
    else:
        z = x.copy()   # start with x's keys and values
        z.update(y)    # modifies z with y's keys and values & returns None
        return z
