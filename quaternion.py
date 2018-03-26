import functools
from math import sqrt

@functools.total_ordering
class Quaternion:
    TOLERANCE = 0.00001
    ROUND = 5

    @classmethod
    def q_normalize(cls, v):
        mag2 = sum(n * n for n in v)
        if abs(mag2 - 1.0) > cls.TOLERANCE:
            mag = sqrt(mag2)
            v = tuple(n / mag for n in v)
        v = tuple(0 if abs(n)<cls.TOLERANCE else round(n,cls.ROUND) for n in v)
        return v

    @classmethod
    def q_conjugate(cls, q):
        w, x, y, z = q
        return Quaternion(w, -x, -y, -z)

    def __init__(self, r=0,i=0,j=0,k=0):
        self.v = self.q_normalize((r,i,j,k))

    def normalize(self):
        mag2 = sum(n * n for n in self.v)
        if abs(mag2 - 1.0) > self.TOLERANCE:
            mag = sqrt(mag2)
            return Quaternion(*[n / mag for n in self.v])
        return self

    def conjugate(self):
        w, x, y, z = self.v
        return Quaternion(w, -x, -y, -z)

    def __format__(self, format_spec):
        return self.__str__()

    def __str__(self):
        st = ''
        if self.v[0]:
            st += str(self.v[0])
            if self.v[1] or self.v[2] or self.v[3]:
                st+='+'
        if self.v[1]:
            st += str(self.v[1])+'i'
            if self.v[2] or self.v[3]:
                st+='+'
        if self.v[2]:
            st += str(self.v[2])+'j'
            if self.v[3]:
                st+='+'
        if self.v[3]:
            st += str(self.v[3])+'k'

        return st if st!='' else '0'

    __repr__ = __str__

    def __hash__(self):
        if not hasattr(self, '_hash'):
            self._hash = hash(self.v)
        return self._hash

    def __mul__(self, other):
        if isinstance(other, Quaternion):
            w1, x1, y1, z1 = self.v
            w2, x2, y2, z2 = other.v
            w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
            x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
            y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
            z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
            return Quaternion(w, x, y, z).normalize()
        return Quaternion(*[(r*other) for r in self.v])

    __rmul__ = __mul__

    def __add__(self, other):
        if isinstance(other, Quaternion):
            return Quaternion(*[self.v[i]+other.v[i] for i in len(self.v)])
        return Quaternion(*[(r+other) for r in self.v])

    __radd__ = __add__

    def __eq__(self, other):
        if other==None:
            return False
        if not isinstance(other, Quaternion):
            return (abs(self.v[1])<self.TOLERANCE and
                    abs(self.v[2])<self.TOLERANCE and
                    abs(self.v[3])<self.TOLERANCE and
                    abs(self.v[0]-other)<self.TOLERANCE)
        return (abs(self.v[0]-other.v[0])<self.TOLERANCE and
                abs(self.v[1]-other.v[1])<self.TOLERANCE and
                abs(self.v[2]-other.v[2])<self.TOLERANCE and
                abs(self.v[3]-other.v[3])<self.TOLERANCE )

    def __lt__(self, other):
        return hash(self) < hash(other)
