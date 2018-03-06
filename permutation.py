import functools

@functools.total_ordering
class Permutation:
    def __init__(self, l, name=None):
        if isinstance(l, dict):
            d = {}
            for i in l:
                if i != l[i]:
                    d[i] = l[i]
            self.d = d
        else:
            self.d = {}
            last = None
            first = None
            for i in l:
                if first==None:
                    first = i
                if last!=None and last != i:
                    self.d[last] = i
                last = i
            if first != None and first != last:
                self.d[last] = first
        if name!=None:
            self.name = name
        else:
            self.name = str(self.d)

    def __format__(self, format_spec):
        if format_spec == '#':
            return self.name
        return str(self)

    def __str__(self):
        return self.name+'='+str(self.d) if self.name!=str(self.d) else str(self.d)

    __repr__ = __str__

    def __hash__(self):
        if not hasattr(self, '_hash'):
            self._hash = hash(frozenset(self.d.items()))
        return self._hash

    def __eq__(self, other):
        if not isinstance(other, Permutation):
            return False
        return self.d == other.d

    def __lt__(self, other):
        if not isinstance(other, Permutation):
            return False
        return hash(self) < hash(other)

    def __getitem__(self, i):
        if i in self.d:
            return self.d[i]
        return i

    def __mul__(self, other):
        d = {}
        if isinstance(other, Permutation):
            for i in other.d:
                nxt = self[other[i]]
                if nxt != i:
                    d[i] = nxt
            for i in self.d:
                if i not in other.d:
                    nxt = self[i]
                    if nxt != i:
                        d[i] = nxt
            return Permutation(d)
        elif other:
            return self
        else:
            return other

    def __add__(self, other):
        if isinstance(other, Permutation):
            return self*other
        elif other:
            return other
        else:
            return self

    __radd__ = __add__
    __rmul__ = __mul__
    __mod__ = __mul__

    def sign(self):
        # Return 1 for even, -1 for odd
        if len(self.d) == 0:
            return 1
        checked = []
        sign = 0
        for i in self.d:
            if i in checked:
                continue
            cur = i
            cycle_length = 1
            while cur not in checked:
                checked.append(cur)
                cur = self.d[cur]
                cycle_length += 1
            sign = (sign+cycle_length)%2
        return -1 if sign else 1
