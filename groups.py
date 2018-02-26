import itertools
import functools
import time
from math import floor
from fractions import gcd

from operations import getop
from matrix import *
from permutation import *

class FiniteGroup:
    def __init__(self, l, op, e=None, name=None):
        self.l = l
        self.op = op
        self.e = e
        self.name = name

        if self.e == None:
            self.identity()

        @functools.lru_cache(len(self.l)**2)
        def _power(i, p):
            term = i
            for k in range(p-1):
                term = self.op(term,i)
            return term
        self.power = _power

        @functools.lru_cache(len(self.l))
        def _order(i):
            for k in range(1, len(self.l)+1):
                if self.power(i,k) == self.e:
                    return k
        self.order = _order

    def identity(self):
        if self.e == None and len(self.l)>0:
            for i in self.l:
                for j in self.l:
                    res = self.op(i,j)
                    if res != j:
                        self.e = None
                        break
                    self.e = i
                if self.e != None:
                    break
        return self.e

    def __hash__(self):
        if not hasattr(self, '_hash'):
            self._hash = hash(frozenset(self.l))
        return self._hash

    def __format__(self, format_spec):
        if format_spec == '#':
            return self.name
        return str(self)

    def __str__(self):
        return '('+', '.join(('['+', '.join([format(i, '#') for i in self.l])+']', str(self.op)))+')'

    __repr__ = __str__

    def __len__(self):
        return len(self.l)

    def __contains__(self, i):
        return i in self.l

    def __eq__(self, other):
        return self.op == other.op and sorted(self.l) == sorted(other.l)

    def __le__(self, other):
        return self.op == other.op and set(self.l).issubset(set(other.l))

    def __lt__(self, other):
        return self<=other and self!=other

    def cayley(self):
        if hasattr(self, '_cayley'):
            return self._cayley
        self._cayley = {}
        for i in self.l:
            for j in self.l:
                self._cayley[(i,j)] = self.op(i,j)
        return self._cayley

    def orders(self):
        if hasattr(self, '_orders'):
            return self._orders
        self._orders = {}
        for i in self.l:
            self._orders[i] = self.order(i)
        return self._orders

    def revorders(self):
        if hasattr(self, '_revorders'):
            return self._revorders
        self._revorders = {}
        for i in self.orders():
            order = self.orders()[i]
            if not order in self._revorders:
                self._revorders[order] = []
            self._revorders[order].append(i)
        return self._revorders

    def cyclic(self):
        if hasattr(self, '_cyclic'):
            return self._cyclic
        c = []
        for i in self.l:
            if sorted(self.l) == sorted([self.power(i,k) for k in range(1,self.order(i)+1)] if self.order(i) else []):
                c.append(i)
        self._cyclic = (len(c) != 0, c)
        return self._cyclic

    def center(self):
        if hasattr(self, '_center'):
            return self._center
        self._center = []
        for i in self.l:
            go = True
            for j in self.l:
                if self.op(i,j) != self.op(j,i):
                    go = False
                    break
            if go:
                self._center.append(i)
        return self._center

    # @functools.lru_cache()
    def centralizer(self, g):
        c = []
        for i in self.l:
            go = True
            for j in g.l:
                if self.op(i,j) != self.op(j,i):
                    go = False
                    break
            if go:
                c.append(i)
        return c

    def subgroup(self, k):
        if len(self.l) == k:
            return self
        if k in self.revorders():
            return GeneratorGroup(self.revorders()[k][0], self.op)
        return NullGroup()

    def subgroups(self, k):
        if len(self.l) == k:
            return [self]
        subs = []
        if k in self.revorders():
            for i in self.revorders()[k]:
                g = GeneratorGroup(i, self.op)
                if not g in subs:
                    subs.append(g)
        return subs if subs!=[] else [NullGroup()]

class NullGroup(FiniteGroup):
    def __init__(self):
        super(NullGroup, self).__init__([], getop('null'), name='NullGroup')

class GeneratorGroup(FiniteGroup):
    def __init__(self, g, op, name=None):
        def power_(i, p):
            term = i
            for k in range(p-1):
                term = op(term,i)
            return term

        if name==None:
            name = '<'+format(g,'#')+'>'
        k = 2
        l = []
        if g != None:
            l = [g]
            p = power_(g,k);
            while not p in l:
                l.append(p)
                k+=1
                p = power_(g,k);
        super(GeneratorGroup, self).__init__(l, op, name=name)

class PermutationGroup(FiniteGroup):
    def __init__(self, l, e=Permutation([]), name=None):
        p = None
        if len(l)>0 and not isinstance(l[0], Permutation):
            p = []
            for i in l:
                p.append(Permutation(i))
        super(PermutationGroup, self).__init__(p if p != None else l, getop('mult'), e, name=name)

class MatrixGroup(FiniteGroup):
    def __init__(self, l, op=getop('mult'), name=None):
        p = None
        if len(l)>0 and not isinstance(l[0], Matrix):
            p = []
            for i in l:
                p.append(Matrix(i))
        e = None
        super(MatrixGroup, self).__init__(p if p!= None else l, op, name=name)

class GL(MatrixGroup):
    def det_case(self, x):
        return x.det()!=0

    def __init__(self, s, g, op=getop('mult'), name=None):
        m = []
        def loop_rec(g, n, l={}, rows=[]):
            if n >= 1:
                if isinstance(g, FiniteGroup):
                    for a in g.l:
                        l[n] = a
                        rows = loop_rec(g, n - 1, l, rows)
                else:
                    for a in g:
                        l[n] = a
                        rows = loop_rec(g, n - 1, l, rows)
            else:
                r = []
                for i in l:
                    r.append(l[i])
                rows.append(r)
            return rows
        rows = loop_rec(g, s)
        matrices = loop_rec(rows, s, {}, [])
        for mat in matrices:
            t = Matrix(mat)
            if self.det_case(t):
                m.append(t)
        super(GL, self).__init__(m, op, name='GL('+str(s)+', '+format(g,'#')+')')

class SL(GL):
    def det_case(self, x):
        return x.det()==1

    def __init__(self, s, g, op=getop('mult')):
        super(SL, self).__init__(s, g, op, name='SL('+str(s)+', '+format(g,'#')+')')

class Aff(MatrixGroup):
    def __init__(self, g, op=getop('mult')):
        m = []
        for a in g.l:
            if a != 0:
                for b in g.l:
                    m.append(Matrix([[a,b],[0,1]]))
        super(Aff, self).__init__(m, op, name='Aff('+format(g,'#')+')')

class S(PermutationGroup):
    def __init__(self, n):
        l = list(itertools.permutations(range(1,n+1)))
        p = []
        for s in l:
            d = {}
            for i in range(n):
                d[i+1] = s[i]
            p.append(d)
        super(S, self).__init__(p, name='S('+str(n)+')')

class D(PermutationGroup):
    def __init__(self, n):
        p = []
        l = []
        r = list(range(1,n+1))
        for i in range(n):
            l.append(r[i:]+r[:i])
        if n > 2:
            r = r[1:]+r[:1]
            for i in range(floor(n/2)):
                t = r[i]
                r[i] = r[n-i-1]
                r[n-i-1] = t
            for i in range(n):
                l.append(list(r[i:])+list(r[:i]))
        name = []
        if n > 0:
            name.append('e')
        if n > 1:
            name.append('f')
        if n > 2:
            name+=['f'+str(i) for i in range(2,n)]+['g', 'gf']+['gf'+str(i) for i in range(2,n)]
        for s in range(len(l)):
            d = {}
            for i in range(n):
                d[i+1] = l[s][i]
            p.append(Permutation(d,name[s]))
        super(D, self).__init__(p, p[0], name='D('+str(n)+')')

class U(FiniteGroup):
    def __init__(self, n):
        super(U, self).__init__([i for i in range(n) if gcd(i,n)==1] if n!=1 else [1], getop('multmod', n), 1, name='U('+str(n)+')')

class Z(FiniteGroup):
    def __init__(self, n):
        super(Z, self).__init__([i for i in range(n)], getop('addmod', n), 0, name='Z('+str(n)+')')

    def subgroup(self, k):
        if len(self.l)%k == 0:
            return GeneratorGroup(int(len(self.l)/k)%len(self.l), self.op)
        return NullGroup()

    def subgroups(self, k):
        return [self.subgroup(k)]

class Zx(FiniteGroup):
    def __init__(self, n):
        super(Zx, self).__init__([i for i in range(1,n)], getop('multmod', n), 1, name='Zx('+str(n)+')')

if __name__ == '__main__':

    def task(group, *args):
        task_start = time.clock()
        g = group(*args)
        print('group:',format(g,'#'),'=',g)
        print('\nidentity:',g.identity(),flush=True)
        if g.identity()!=None:
            print('\ncyclic:',g.cyclic())
            print('\norders:',g.orders())
            print('\nsubgroups:')
            for i in range(1,len(g)+1):
                subs = g.subgroups(i)
                if subs[0] != NullGroup():
                    print(i)
                    for s in subs:
                        print('\t:',s)
        task_end = time.clock()
        print('\ntask_time =',task_end - task_start,'s\n')
        print('*****************************************\n',flush=True)

    start = time.clock()

    # Main code starts here

    task(SL, 3, Z(2), getop('matrixmod', 2))
    task(Aff, D(2))
    task(D, 5)
    task(S, 5)
    task(U, 283)

    # Main code ends here

    end = time.clock()
    print('\ntotal_time =',end - start,'s')
