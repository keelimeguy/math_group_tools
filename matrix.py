import functools
from numpy import matrix

@functools.total_ordering
class Matrix:
    def __init__(self, m):
        if isinstance(m, matrix):
            self.m = m
        else:
            self.m = matrix(m)

    def rows(self):
        return self.m.shape[0]

    def cols(self):
        return self.m.shape[1]

    def cut(self, r, c):
        l = self.m.tolist()
        m = []
        for i in range(self.rows()):
            if i != r:
                n = []
                for j in range(self.cols()):
                    if j != c:
                        n.append(l[i][j])
                m.append(n)
        return Matrix(m)

    def det(self):
        if hasattr(self, '_det'):
            return self._det
        if self.rows()!=self.cols():
            return None
        if self.rows()==1:
            return self.m.tolist()[0][0]
        self._det = 0
        for c in range(self.cols()):
            self._det += ((-1)**(c%2)) * self.m.tolist()[0][c] * self.cut(0, c).det()
        return self._det

    def __format__(self, format_spec):
        return '['+', '.join(['['+', '.join([format(c, format_spec) for c in r])+']' for r in self.m.tolist()])+']'

    def __str__(self):
        return str(self.m.tolist())

    __repr__ = __str__

    def __hash__(self):
        if not hasattr(self, '_hash'):
            l = self.m.tolist()
            s = []
            for i in l:
                for e in range(len(i)):
                    s.append((e+1)*hash(frozenset(i))*hash(i[e]))
            self._hash = hash(frozenset(s))
        return self._hash

    def __mul__(self, other):
        if isinstance(other, Matrix):
            return Matrix(self.m*other.m)
        return Matrix([[c*other for c in r] for r in self.m.tolist()])

    __rmul__ = __mul__

    def __add__(self, other):
        if isinstance(other, Matrix):
            return Matrix(self.m+other.m)
        return Matrix([[c+other for c in r] for r in self.m.tolist()])

    __radd__ = __add__

    def __eq__(self, other):
        if not isinstance(other, Matrix):
            return False
        if self.rows()!=other.rows() or self.cols()!=other.cols():
            return False
        s_m = self.m.tolist()
        o_m = other.m.tolist()
        for r in range(self.rows()):
            for c in range(self.cols()):
                if s_m[r][c] != o_m[r][c]:
                    return False
        return True

    def __lt__(self, other):
        return hash(self) < hash(other)
