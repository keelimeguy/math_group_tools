import functools

from matrix import *
from permutation import *

class Operation:
    def __init__(self, op, name):
        self.name = name
        self.op = op

    def __str__(self):
        return self.name

    __repr__ = __str__

    def __call__(self, *args):
        return self.op(*args)

ops = { 'add':Operation(lambda x,y: x+y, 'add'),
        'mult':Operation(lambda x,y: x*y, 'mult'),
        'matrixmod':(lambda mod: Operation(lambda x,y: Matrix([[c%mod for c in r] for r in (x*y).m.tolist()]), 'matrixmod'+str(mod)), {}),
        'matrixelement':(lambda eop: Operation(lambda x,y: Matrix([[eop(x.m.tolist()[r][c],y.m.tolist()[r][c]) for c in range(len(x.m.tolist()[r]))] for r in range(len(x.m.tolist()))]), 'matrixelement'+str(eop)), {}),
        'addmod':(lambda mod: Operation(lambda x,y: (x+y)%mod, 'addmod'+str(mod)), {}),
        'multmod':(lambda mod: Operation(lambda x,y: (x*y)%mod, 'multmod'+str(mod)), {})}

def getop(type, n=None):
    if type in ops:
        if n!=None and not isinstance(ops[type], Operation):
            if isinstance(n, Operation):
                if not str(n) in ops[type][1]:
                    ops[type][1][str(n)] = ops[type][0](n)
                return ops[type][1][str(n)]
            else:
                if not n in ops[type][1]:
                    ops[type][1][n] = ops[type][0](n)
                return ops[type][1][n]
        else:
            return ops[type]

if __name__ == '__main__':
    o = getop('matrixmod', 7)
    print(o, o(Matrix([[1,2],[3,4]]), Matrix([[1,2],[3,4]])))
    o = getop('null')
    print(o)
