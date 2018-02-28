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
        return self.op == other.op and self.sorted() == other.sorted()

    def __le__(self, other):
        return self.op == other.op and set(self.l).issubset(set(other.l))

    def __lt__(self, other):
        return self<=other and self!=other

    def sorted(self):
        if not hasattr(self, '_sorted'):
            self._sorted = sorted(self.l)
        return list(self._sorted)

    def cayley(self):
        if not hasattr(self, '_cayley'):
            self._cayley = {}
            l = self.sorted()
            r = self.sorted()
            for i in l:
                for j in r:
                    # print(i,j,flush=True)
                    res = self.op(i,j)
                    self._cayley[(i,j)] = res
                    if hasattr(self, '_abelian') and self._abelian:
                        self._cayley[(j,i)] = res
                if hasattr(self, '_abelian') and self._abelian:
                    r.remove(i)
        return self._cayley

    def orders(self):
        if not hasattr(self, '_orders'):
            self._orders = {}
            for i in self.l:
                self._orders[i] = self.order(i)
        return self._orders

    def revorders(self):
        if not hasattr(self, '_revorders'):
            self._revorders = {}
            for i in self.orders():
                order = self.orders()[i]
                if not order in self._revorders:
                    self._revorders[order] = []
                self._revorders[order].append(i)
        return self._revorders

    def cyclic(self):
        if not hasattr(self, '_cyclic'):
            c = []
            for i in self.l:
                # print(i,flush=True)
                if self.sorted() == sorted([self.power(i,k) for k in range(1,self.order(i)+1)] if self.order(i) else []):
                    c.append(i)
            self._cyclic = (len(c) != 0, c)
            if self._cyclic[0] and not hasattr(self, '_abelian'):
                self._abelian = True
        return self._cyclic

    def abelian(self):
        if not hasattr(self, '_abelian'):
            self._abelian = True
            for i in self.l:
                for j in self.l:
                    # print(i,j,flush=True)
                    if self.op(i,j) != self.op(j,i):
                        self._abelian = False
                        break
                if not self._abelian:
                    break
        return self._abelian

    def center(self):
        if not hasattr(self, '_center'):
            if hasattr(self, '_abelian') and self._abelian:
                self._center = self.sorted()
            else:
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
        if hasattr(self, '_cyclic') and self._cyclic[0]:
            if k in self.revorders():
                subs.append(GeneratorGroup(self.revorders()[k][0], self.op))
        else:
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

        self._abelian = True
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
        super(PermutationGroup, self).__init__(p if p != None else l, getop('mult', cache=len(l)**2), e, name=name)

class MatrixGroup(FiniteGroup):
    def __init__(self, l, op=getop('mult'), name=None):
        p = None
        if len(l)>0 and not isinstance(l[0], Matrix):
            p = []
            for i in l:
                p.append(Matrix(i))
        e = None
        super(MatrixGroup, self).__init__(p if p!= None else l, op, name=name)

class M(MatrixGroup):
    def det_case(self, x):
        return True

    def __init__(self, s, g, op=None, name=None):
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
        if name == None:
            name = 'M('+str(s)+', '+format(g,'#')+')'
        if op==None:
            op = getop('mult', cache=len(m)**2)
        super(M, self).__init__(m, op, name)

class GL(M):
    def det_case(self, x):
        return x.det()!=0

    def __init__(self, s, g, op=None):
        super(GL, self).__init__(s, g, op, name='GL('+str(s)+', '+format(g,'#')+')')

class SL(M):
    def det_case(self, x):
        return x.det()==1

    def __init__(self, s, g, op=None):
        super(SL, self).__init__(s, g, op, name='SL('+str(s)+', '+format(g,'#')+')')

class Aff(MatrixGroup):
    def __init__(self, g, op=None):
        m = []
        for a in g.l:
            if a != 0:
                for b in g.l:
                    m.append(Matrix([[a,b],[0,1]]))
        if op==None:
            op = getop('mult', (len(m)**2))
        super(Aff, self).__init__(m, op, name='Aff('+format(g,'#')+')')

class S(PermutationGroup):
    def __init__(self, n):
        if n>0 and n<=2:
            self._abelian = True
        else:
            self._abelian = False
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
            self._abelian = True
            name.append('e')
        if n > 1:
            name.append('f')
        if n > 2:
            self._abelian = False
            name+=['f'+str(i) for i in range(2,n)]+['g', 'gf']+['gf'+str(i) for i in range(2,n)]
        for s in range(len(l)):
            d = {}
            for i in range(n):
                d[i+1] = l[s][i]
            p.append(Permutation(d,name[s]))
        super(D, self).__init__(p, p[0], name='D('+str(n)+')')

class U(FiniteGroup):
    def __init__(self, n):
        self._abelian = True
        super(U, self).__init__([i for i in range(n) if gcd(i,n)==1] if n!=1 else [1], getop('multmod', n, cache=0), 1 if n>0 else None, name='U('+str(n)+')')

class Z(FiniteGroup):
    def __init__(self, n):
        self._abelian = True
        super(Z, self).__init__([i for i in range(n)], getop('addmod', n, cache=0), 0 if n>0 else None, name='Z('+str(n)+')')

    def subgroup(self, k):
        if len(self.l)%k == 0:
            return GeneratorGroup(int(len(self.l)/k)%len(self.l), self.op)
        return NullGroup()

    def subgroups(self, k):
        return [self.subgroup(k)]

class Zx(FiniteGroup):
    def __init__(self, n):
        self._abelian = True
        super(Zx, self).__init__([i for i in range(1,n)], getop('multmod', n, cache=0), 1 if n>1 else None, name='Zx('+str(n)+')')

if __name__ == '__main__':
    def task(group, *args):
        task_start = time.clock()
        g = group(*args)
        print('group:',format(g,'#'),'=',g)
        print('\nlength:',len(g.l))
        print('\nidentity:',g.identity(),flush=True)
        if g.identity()!=None:
            print('\ncyclic:',g.cyclic(),flush=True)
            print('\norders:',g.orders(),flush=True)
            print('\nabelian:',g.abelian(),flush=True)
            print('\ncenter:',g.center(),flush=True)
            print('\ncayley:',g.cayley(),flush=True)
            print('\nsubgroups:')
            for i in range(1,len(g)+1):
                subs = g.subgroups(i)
                if subs[0] != NullGroup():
                    print(i)
                    for s in subs:
                        print('\t:',s,flush=True)
            print('\ncache:')
            print('\tpower:',g.power.cache_info())
            print('\torder:',g.power.cache_info())
            print('\top:',g.op.cache_info())
        task_end = time.clock()
        print('\ntask_time =',task_end - task_start,'s\n')
        print('*****************************************\n',flush=True)

    start = time.clock()

    # Main code starts here
    n=2
    g=4
    z = Z(g)
    task(M, n, z, getop('matrixelement', z.op, id=1, cache=(g**(n**2))**2))
    task(GL, n, z, getop('matrixelement', z.op, id=1))
    task(SL, n, z, getop('matrixelement', z.op, id=1))
    n=2
    g=3
    z = Z(g)
    task(M, n, z, getop('matrixmod', g, id=1, cache=(g**(n**2))**2))
    task(GL, n, z, getop('matrixmod', g, id=1))
    task(SL, n, z, getop('matrixmod', g, id=1))
    task(Aff, D(6))
    task(D, 5)
    task(S, 5)
    task(U, 263)

    # Main code ends here

    end = time.clock()
    print('\ntotal_time =',end - start,'s')
