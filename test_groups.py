import unittest
from groups import *
from operations import getop
from matrix import *
from permutation import *

class GroupsTestCase(unittest.TestCase):
    def test_matrixgroups(self):
        a = Aff(Z(4),getop('matrixmod', 4))
        self.assertEqual(a.orders()[Matrix([[3,3],[0,1]])], 2)
        self.assertEqual(a.orders()[Matrix([[2,1],[0,1]])], None)
        z = Z(3)
        s = SL(3, z, getop('matrixelement', z.op))
        g = GL(3, z, getop('matrixelement', z.op))
        self.assertEqual(s<g, True)
        self.assertEqual(s<=g, True)
        self.assertEqual(s>g, False)
        self.assertEqual(g<s, False)
        self.assertEqual(g>s, True)
        self.assertEqual(g>=s, True)

    def test_D(self):
        self.assertEqual(D(4), PermutationGroup([
            Permutation([]),
            Permutation([1, 2, 3, 4]),
            Permutation([1, 3])*Permutation([2, 4]),
            Permutation([1, 4, 3, 2]),
            Permutation([2, 4]),
            Permutation([1, 4])*Permutation([2, 3]),
            Permutation([1, 3]),
            Permutation([1, 2])*Permutation([3, 4])]))
        self.assertEqual(D(4)<S(4), True)
        self.assertEqual(D(3), S(3))

    def test_permutations(self):
        self.assertEqual(Permutation([1, 3, 4])[3], 4)
        self.assertEqual(Permutation([1, 3, 4])[4], 1)
        self.assertEqual(Permutation([1, 3]),
            Permutation([3, 1]))
        self.assertEqual(Permutation([1, 3]),
            Permutation([2, 1])*
            Permutation([2, 3])*
            Permutation([2, 1]))

    def test_generator(self):
        for i in range(5,20):
            self.assertEqual(Z(i),
                GeneratorGroup(1 if i>1 else (0 if i==1 else None), getop('addmod', i)))

    def test_combinational_matrix(self):
        self.assertEqual(Aff(D(2)), GeneratorGroup(Matrix([
            [Permutation([1,2]),Permutation([1,2])],
            [0,1]]), Aff(D(2)).op))
        self.assertEqual(Aff(D(2)), GeneratorGroup(Matrix([
            [Permutation([1,2]),Permutation([])],
            [0,1]]), Aff(D(2)).op))

if __name__ == '__main__':
    # About 15 seconds, due to test_matrixgroups
    unittest.main()
