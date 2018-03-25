import itertools
import functools
import argparse
import time
import re
from math import floor
from fractions import gcd

from operations import getop
from matrix import *
from permutation import *

class FiniteGroup:
    def print_help():
        print('FiniteGroup arguments:')
        print('\t','l = list of elements')
        print('\t','op = group operation')
        print('\t','?e = identity element')
        print('\t','?name = name of the group')

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
        if format_spec == '#' and self.name != None:
            return self.name
        return str(self)

    def __str__(self):
        return '('+', '.join(('['+', '.join([format(i, '#') for i in self.l])+']', str(self.op)))+')'

    __repr__ = __str__

    def __len__(self):
        return len(self.l)

    def __iter__(self):
        for e in self.l:
            yield e

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
                if sorted(self._center) == self.sorted():
                    self._abelian = True
                else:
                    self._abelian = False
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

    def lcosets(self, h):
        if not h<=self:
            return 'Warning: '+format(h,'#')+' is not a subset of '+format(self,'#')
        cosets = [h.l]
        remaining = []
        for i in self.sorted():
            if i not in cosets[0]:
                remaining.append(i)
        while remaining:
            i+=1
            a = remaining.pop()
            cosets.append([])
            for n in cosets[0]:
                o = self.op(a,n)
                for j, r in enumerate(remaining):
                    if r == o:
                        del remaining[j]
                        break
                cosets[i].append(o)
        return cosets

    def rcosets(self, h):
        if not h<=self:
            return 'Warning: '+h.name+' is not a subset of '+self.name
        cosets = [h.l]
        remaining = []
        for i in self.sorted():
            if i not in cosets[0]:
                remaining.append(i)
        i = 0
        while remaining:
            i+=1
            a = remaining.pop()
            cosets.append([])
            for n in cosets[0]:
                o = self.op(n,a)
                for j, r in enumerate(remaining):
                    if r == o:
                        del remaining[j]
                        break
                cosets[i].append(o)
        return cosets

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
        return subs if subs else [NullGroup()]

class NullGroup(FiniteGroup):
    def print_help():
        print('NullGroup arguments: None')

    def __init__(self):
        super(NullGroup, self).__init__([], getop('null'), name='NullGroup')

class GeneratorGroup(FiniteGroup):
    def print_help():
        print('GeneratorGroup arguments:')
        print('\t','g = generator element or list of generator elements')
        print('\t','op = group operation')
        print('\t','?name = name of the group')

    def __init__(self, g, op, name=None):
        def power_(i, p):
            term = i
            for k in range(p-1):
                term = op(term,i)
            return term

        d = {}
        l_g = [g]
        if isinstance(g,list):
            l_g = g
            if name==None:
                name = '<'+','.join([format(i,'#') for i in g])+'>'
        elif name==None:
            self._abelian = True
            name = '<'+format(g,'#')+'>'
        for g in l_g:
            k = 2
            if g not in d:
                d[g] = [g]
                p = power_(g,k);
                while not p in d[g]:
                    d[g].append(p)
                    k+=1
                    p = power_(g,k);

        l = []
        for k in d:
            for i in d[k]:
                if not i in l:
                    l.append(i)
                for j in d:
                    o = op(i,j)
                    if not o in l:
                        l.append(o)

        super(GeneratorGroup, self).__init__(l, op, name=name)

class MatrixGeneratorGroup(GeneratorGroup):
    def print_help():
        print('MatrixGeneratorGroup arguments:')
        print('\t','m = generator element, a matrix or list of matrices')
        print('\t','op = group operation')
        print('\t','?name = name of the group')

    def __init__(self, m, op, name=None):
        g = m
        if not isinstance(m, Matrix):
            if len(m[0])>0 and isinstance(m[0][0], list):
                g = []
                for h in m:
                    g.append(Matrix(h))
            else:
                g = Matrix(m)
        super(MatrixGeneratorGroup, self).__init__(g, op, name=name)

class PermutationGeneratorGroup(GeneratorGroup):
    def print_help():
        print('PermutationGeneratorGroup arguments:')
        print('\t','p = generator element, a permutation')
        print('\t','?pname = name for the generator element, when not already named')
        print('\t','?name = name of the group')

    def __init__(self, p, pname=None, name=None):
        g = p
        if not isinstance(p, Permutation):
            if len(p[0])>0 and isinstance(p[0][0], list):
                g = []
                for h in p:
                    g.append(Permutation(h))
            else:
                g = Permutation(p, pname)
        super(PermutationGeneratorGroup, self).__init__(g, getop('mult'), name=name)

class PermutationGroup(FiniteGroup):
    def print_help():
        print('PermutationGroup arguments:')
        print('\t','l = list of elements')
        print('\t','?e = identity element')
        print('\t','?name = name of the group')

    def __init__(self, l, e=Permutation([]), name=None):
        p = None
        if len(l)>0 and not isinstance(l[0], Permutation):
            p = []
            for i in l:
                p.append(Permutation(i))
        super(PermutationGroup, self).__init__(p if p != None else l, getop('mult', cache=len(l)**2), e, name=name)

class MatrixGroup(FiniteGroup):
    def print_help():
        print('MatrixGroup arguments:')
        print('\t','l = list of elements')
        print('\t','?op = group operation')
        print('\t','?name = name of the group')

    def __init__(self, l, op=getop('mult'), name=None):
        p = None
        if len(l)>0 and not isinstance(l[0], Matrix):
            p = []
            for i in l:
                p.append(Matrix(i))
        e = None
        super(MatrixGroup, self).__init__(p if p!= None else l, op, name=name)

class M(MatrixGroup):
    def print_help():
        print('M arguments:')
        print('\t','s = square size of matrices')
        print('\t','g = group or list of elements')
        print('\t','?op = group operation')
        print('\t','?name = name of the group')

    def det_case(self, x):
        return True

    def __init__(self, s, g, op=None, name=None):
        s = int(s)
        m = []
        def loop_rec(g, n, l={}, rows=[]):
            if n >= 1:
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
            name = 'M('+str(s)+', '+(format(g,'#') if isinstance(g, FiniteGroup) else str(g))+')'
        if op==None:
            op = getop('mult', cache=len(m)**2)
        super(M, self).__init__(m, op, name)

class GL(M):
    def print_help():
        print('GL arguments:')
        print('\t','s = square size of matrices')
        print('\t','g = group or list of elements')
        print('\t','?op = group operation')

    def det_case(self, x):
        return x.det()!=0

    def __init__(self, s, g, op=None):
        s = int(s)
        super(GL, self).__init__(s, g, op, name='GL('+str(s)+', '+(format(g,'#') if isinstance(g, FiniteGroup) else str(g))+')')

class SL(M):
    def print_help():
        print('SL arguments:')
        print('\t','s = square size of matrices')
        print('\t','g = group or list of elements')
        print('\t','?op = group operation')

    def det_case(self, x):
        return x.det()==1

    def __init__(self, s, g, op=None):
        s = int(s)
        super(SL, self).__init__(s, g, op, name='SL('+str(s)+', '+(format(g,'#') if isinstance(g, FiniteGroup) else str(g))+')')

class Aff(MatrixGroup):
    def print_help():
        print('Aff arguments:')
        print('\t','g = group or list of elements')
        print('\t','?op = group operation')

    def __init__(self, g, op=None):
        m = []
        for a in g:
            if a != 0:
                for b in g:
                    m.append(Matrix([[a,b],[0,1]]))
        if op==None:
            op = getop('mult', (len(m)**2))
        super(Aff, self).__init__(m, op, name='Aff('+(format(g,'#') if isinstance(g, FiniteGroup) else str(g))+')')

class S(PermutationGroup):
    def print_help():
        print('S arguments:')
        print('\t','n = group number')

    def __init__(self, n):
        n = int(n)
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

class A(PermutationGroup):
    def print_help():
        print('A arguments:')
        print('\t','n = group number')

    def __init__(self, n):
        n = int(n)
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
            perm = Permutation(d)
            if perm.sign()>0:
                p.append(perm)
        super(A, self).__init__(p, name='A('+str(n)+')')

class D(PermutationGroup):
    def print_help():
        print('D arguments:')
        print('\t','n = group number')

    def __init__(self, n):
        n = int(n)
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
    def print_help():
        print('U arguments:')
        print('\t','n = group number')

    def __init__(self, n):
        n = int(n)
        self._abelian = True
        super(U, self).__init__([i for i in range(n) if gcd(i,n)==1] if n!=1 else [1], getop('multmod', n, cache=0), 1 if n>0 else None, name='U('+str(n)+')')

class Z(FiniteGroup):
    def print_help():
        print('Z arguments:')
        print('\t','n = group number')

    def __init__(self, n):
        n = int(n)
        self._abelian = True
        super(Z, self).__init__([i for i in range(n)], getop('addmod', n, cache=0), 0 if n>0 else None, name='Z('+str(n)+')')

    def subgroup(self, k):
        if len(self.l)%k == 0:
            return GeneratorGroup(int(len(self.l)/k)%len(self.l), self.op)
        return NullGroup()

    def subgroups(self, k):
        return [self.subgroup(k)]

class Zx(FiniteGroup):
    def print_help():
        print('Zx arguments:')
        print('\t','n = group number')

    def __init__(self, n):
        n = int(n)
        self._abelian = True
        super(Zx, self).__init__([i for i in range(1,n)], getop('multmod', n, cache=0), 1 if n>1 else None, name='Zx('+str(n)+')')





if __name__ == '__main__':
    print()

    def append_required_length(nmin,nmax):
        class AppendRequiredLength(argparse.Action):
            def __call__(self, parser, args, values, option_string=None):
                if not nmin<=len(values)<=nmax:
                    msg='Argument "{f}" requires between {nmin} and {nmax} arguments'.format(
                        f=self.dest,nmin=nmin,nmax=nmax)
                    raise argparse.ArgumentTypeError(msg)
                if hasattr(args, self.dest):
                    current_values = getattr(args, self.dest)
                    try:
                        current_values.extend([values])
                    except AttributeError:
                        current_values = [values]
                    finally:
                        setattr(args, self.dest, current_values)
                else:
                    setattr(args, self.dest, [values])
        return AppendRequiredLength

    parser = argparse.ArgumentParser(
        description='Find properties/desciptions of a given group. Sample input:\
                    \'python groups.py -g S 3 -t abelian orders "centralizer{S,3}" "lcosets{PermutationGroup,[[1,2]]}" "rcosets{PermutationGroup,[[1,2]]}" cache\'')
    parser.add_argument('-g', '--group', action=append_required_length(1,2), nargs='+', metavar=('group','args'),
                        help='The type of group to use, then comma separated arguments to the group.')
    parser.add_argument('-i', '--info', action="store_true",
                        help="Give info on available groups or arguments to a group.\
                            e.g. \'-i -g <group>\'")
    parser.add_argument('-t', '--task', action="append", nargs="+",
                        help='The subtasks to perform, as space separated list.\
                            Default is \'-t cyclic orders abelian center cayley subgroups cache\'.\
                            Possible tasks: cyclic, orders, abelian, center, centralizer{group}, lcosets{group}, rcosets{group}, cayley, subgroups, cache')
    parser.add_argument('-n', '--notask', action="append", nargs="+",
                        help='The subtasks to skip, as space separated list. Ignored if --task is used.\
                            Default is \'-t cyclic orders abelian center cayley subgroups cache\', the --notask argument\
                            allows you to skip any of these default tasks.\
                            Possible tasks: cyclic, orders, abelian, center, cayley, subgroups, cache')

    args = parser.parse_args()

    group_type_list = {'Zx':Zx, 'Z':Z, 'U':U, 'A':A, 'D':D, 'S':S, 'Aff':Aff, 'SL':SL, 'GL':GL, 'M':M,
                        'FiniteGroup':FiniteGroup, 'PermutationGeneratorGroup':PermutationGeneratorGroup,
                        'MatrixGeneratorGroup':MatrixGeneratorGroup, 'GeneratorGroup':GeneratorGroup,
                        'PermutationGroup':PermutationGroup, 'MatrixGroup':MatrixGroup}

    if args.info:
        if args.group:
            print('Arguments are comma separated')
            print('\tlists are given as, comma separated within []')
            print('\tgroups are given as, type then arguments comma separated within {}')
            print('\toperations are given as, type then arguments comma separated within <>')
            print('\te.g. \'python groups.py -g S 4 -g M \"2,{Z,2},<matrixelement,<addmod,2>,cache=256>\"\'')
            for g in args.group:
                if g[0] in group_type_list:
                    group_type_list[g[0]].print_help()
                else:
                    print('Warning: Group type', g[0], 'not available.')
        else:
            print("Available groups:", sorted(list(group_type_list)))
        exit(0)

    if args.group == None:
        print('Try running with -h next time.')
        exit(0)

    def subgroup_task(g):
        print('\nsubgroups:')
        for i in range(1,len(g)+1):
            subs = g.subgroups(i)
            if subs[0] != NullGroup():
                print(i)
                for s in subs:
                    print('\t:',s,flush=True)

    def cache_task(g):
        print('\ncache:')
        print('\tpower:',g.power.cache_info())
        print('\torder:',g.power.cache_info())
        print('\top:',g.op.cache_info())

    task_list = ['cyclic', 'orders', 'abelian', 'center', 'centralizer', 'lcosets', 'rcosets', 'cayley', 'subgroups', 'cache']
    default_task_list = ['cyclic', 'orders', 'abelian', 'center', 'cayley', 'subgroups', 'cache']
    task_dict = {
        'cyclic':(lambda g: print('\ncyclic:',g.cyclic(),flush=True)),
        'orders':(lambda g: print('\norders:',g.orders(),flush=True)),
        'abelian':(lambda g: print('\nabelian:',g.abelian(),flush=True)),
        'center':(lambda g: print('\ncenter:',g.center(),flush=True)),
        'centralizer':(lambda g, h: print('\ncentralizer('+format(h,'#')+'):',g.centralizer(h),flush=True)),
        'lcosets':(lambda g, h: print('\nlcosets('+format(h,'#')+'):',g.lcosets(h),flush=True)),
        'rcosets':(lambda g, h: print('\nrcosets('+format(h,'#')+'):',g.rcosets(h),flush=True)),
        'cayley':(lambda g: print('\ncayley:',g.cayley(),flush=True)),
        'subgroups':subgroup_task,
        'cache':cache_task,
    }

    tasks_to_perform = []
    if args.task:
        for t in args.task[0]:
            t = re.compile("{|}").split(t)
            if t[0] in task_list:
                tasks_to_perform.append(t)
    elif args.notask:
        for t in default_task_list:
            if t not in args.notask[0]:
                tasks_to_perform.append([t])
    else:
        for t in default_task_list:
            tasks_to_perform.append([t])

    def task(group, *args, **kwargs):
        task_start = time.clock()
        g = group(*args, **kwargs)
        print('group:',format(g,'#'),'=',g)
        print('\nlength:',len(g))
        print('\nidentity:',g.identity(),flush=True)
        if g.identity()!=None:
            for task in tasks_to_perform:
                if len(task)>1:
                    if ',' in task[1]:
                        h = parse_group([''.join(task[1][:task[1].index(',')]), ''.join(task[1][task[1].index(',')+1:])])
                        task_dict[task[0]](g, h[0](*h[1], **h[2]))
                    else:
                        h = parse_group([task[1],""])
                        task_dict[task[0]](g, h[0](*h[1], **h[2]))
                else:
                    task_dict[task[0]](g)
        task_end = time.clock()
        print('\ntask_time =',task_end - task_start,'s\n')
        print('*****************************************\n',flush=True)
        return task_end - task_start

    def parse_group(g):
        if g[0] not in group_type_list:
            print('Warning: Group type', g[0], 'not available.')
            return None
        else:
            skip = False
            in_equality = []
            scope=[]
            group_args = []
            if len(g) > 1:
                g_args = re.compile("(<|>|{|}|\[|\]|,|=)").split(g[1])
                for i, arg in enumerate(g_args):
                    if skip or arg in [',', '']:
                        skip = False
                        continue
                    elif arg in ['<', '{', '[']:
                        scope = [(arg,[])] + scope
                    elif scope:
                        if scope[0][0] == '{':
                            if arg == '}':
                                new_scope_arg = [[],{}]
                                for o in scope[0][1]:
                                    if isinstance(o,dict):
                                        for k in o:
                                            new_scope_arg[1][k] = o[k]
                                    else:
                                        new_scope_arg[0].append(o)
                                if len(scope)>1:
                                    scope[1][1].append(new_scope_arg[0][0](*new_scope_arg[0][1:],**new_scope_arg[1]))
                                else:
                                    if in_equality:
                                        group_args.append({in_equality[0]: new_scope_arg[0][0](*new_scope_arg[0][1:],**new_scope_arg[1])})
                                    else:
                                        group_args.append(new_scope_arg[0][0](*new_scope_arg[0][1:],**new_scope_arg[1]))
                                scope.pop(0)
                            elif arg == '=':
                                scope[0][1].pop()
                                scope[0][1].append({g_args[i-1]: g_args[i+1]})
                                skip = True
                            elif scope[0][1]:
                                scope[0][1].append(arg)
                            elif arg in group_type_list:
                                scope[0][1].append(group_type_list[arg])
                            else:
                                print('Warning: Group type', arg, 'not available.')
                                return None
                        elif scope[0][0] == '[':
                            if arg == ']':
                                if len(scope)>1:
                                    scope[1][1].append(scope[0][1])
                                else:
                                    if in_equality:
                                        group_args.append({in_equality[0]: scope[0][1]})
                                    else:
                                        group_args.append(scope[0][1])
                                scope.pop(0)
                            else:
                                # Parse list elements as ints until better solution found
                                scope[0][1].append(int(arg))
                        elif scope[0][0] == '<':
                            if arg == '>':
                                new_scope_arg = [[],{}]
                                for o in scope[0][1]:
                                    if isinstance(o,dict):
                                        for k in o:
                                            new_scope_arg[1][k] = o[k]
                                    else:
                                        new_scope_arg[0].append(o)
                                if len(scope)>1:
                                    scope[1][1].append(getop(new_scope_arg[0][0],*new_scope_arg[0][1:],**new_scope_arg[1]))
                                else:
                                    if in_equality:
                                        group_args.append({in_equality[0]: getop(new_scope_arg[0][0],*new_scope_arg[0][1:],**new_scope_arg[1])})
                                    else:
                                        group_args.append(getop(new_scope_arg[0][0],*new_scope_arg[0][1:],**new_scope_arg[1]))
                                scope.pop(0)
                            elif arg == '=':
                                scope[0][1].pop()
                                scope[0][1].append({g_args[i-1]: g_args[i+1]})
                                skip = True
                            else:
                                scope[0][1].append(arg)
                    elif arg == '=':
                        group_args.pop()
                        if g_args[i+1] not in ['<', '{', '[', '']:
                            group_args.append({g_args[i-1]: g_args[i+1]})
                            skip = True
                        else:
                            in_equality = [g_args[i-1]]
                    elif arg:
                        group_args.append(arg)
                    if i == len(g_args)-1 and (scope or skip):
                        return None
            new_group_arg = [[],{}]
            for o in group_args:
                if isinstance(o,dict):
                    for k in o:
                        new_group_arg[1][k] = o[k]
                else:
                    new_group_arg[0].append(o)
            return (group_type_list[g[0]], new_group_arg[0], new_group_arg[1])

    groups = []
    for g in args.group:
        new_group = parse_group(g)
        if new_group:
            groups.append(new_group)

    total_task_time = 0
    for g in groups:
        total_task_time += task(g[0], *g[1], **g[2])

    print('\ntotal_time =',total_task_time,'s')
