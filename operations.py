import functools

from matrix import *
from permutation import *

class Operation:
    def __init__(self, op, name, cache=128):
        self.op = op
        self.name = name
        self.cache = cache
        if cache!=0:
            @functools.lru_cache(cache)
            def _call(*args):
                return self.op(*args)
        else:
            def _call(*args):
                return self.op(*args)
        self.call = _call

    def __str__(self):
        return self.name

    __repr__ = __str__

    def __eq__(self, other):
        if not isinstance(other, Operation):
            return False
        return self.name == other.name

    def __call__(self, *args):
        return self.call(*args)

    def cache_info(self):
        if self.cache!=0:
            return self.call.cache_info()
        else:
            return 'None'

ops = { 'add':(lambda *args: Operation(lambda x,y: x+y, 'add', *args),{}),
        'mult':(lambda *args: Operation(lambda x,y: x*y, 'mult', *args),{}),
        'matrixmod':(lambda mod,*args: Operation(lambda x,y: Matrix([[c%mod for c in r] for r in (x*y).m.tolist()]), 'matrixmod'+str(mod), *args),{}),
        'matrixelement':(lambda eop,*args: Operation(lambda x,y: Matrix([[eop(x.m.tolist()[r][c],y.m.tolist()[r][c]) for c in range(len(x.m.tolist()[r]))] for r in range(len(x.m.tolist()))]), 'matrixelement'+str(eop), *args),{}),
        'addmod':(lambda mod,*args: Operation(lambda x,y: (x+y)%mod, 'addmod'+str(mod), *args),{}),
        'multmod':(lambda mod,*args: Operation(lambda x,y: (x*y)%mod, 'multmod'+str(mod), *args),{})}

def getop(type, *args, **kwargs):
    op = None
    if type in ops:
        if 'id' in kwargs:
            if kwargs['id'] in ops[type][1]:
                return ops[type][1][kwargs['id']]
        if 'cache' in kwargs:
            op = ops[type][0](*args, kwargs['cache'])
        else:
            op = ops[type][0](*args)
        if 'id' in kwargs:
            ops[type][1][kwargs['id']] = op
    return op
