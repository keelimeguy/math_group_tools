import unittest
from groups import *
from operations import *
from matrix import *
from permutation import *
from quaternion import *

class GroupsTestCase(unittest.TestCase):
    def test_operations(self):
        o = getop('mult')
        self.assertEqual(o, getop('mult'))
        o.name = 'x'
        self.assertEqual(getop('mult').name, 'mult')
        o = getop('mult', id='test_operations')
        self.assertEqual(o, getop('mult'))
        o.name = 'x'
        self.assertEqual(getop('mult').name, 'mult')
        self.assertEqual(getop('mult', id='test_operations').name, 'x')
        self.assertEqual(getop('matrixmod', 7)(Matrix([[1,2],[3,4]]), Matrix([[1,2],[3,4]])), Matrix([[0,3],[1,1]]))

    def test_consistency(self):
        g = Z(3)
        l1 = g.sorted()
        l2 = g.sorted()
        self.assertEqual(l1, l2)
        l2.remove(g.e)
        self.assertEqual(l1, g.sorted())
        self.assertEqual(l2 == g.sorted(), False)

    def test_matrixgroups(self):
        a = Aff(Z(4),getop('matrixmod', 4))
        self.assertEqual(a.orders()[Matrix([[3,3],[0,1]])], 2)
        self.assertEqual(a.orders()[Matrix([[2,1],[0,1]])], None)

    def test_comparisons(self):
        z = Z(3)
        m = M(2, z, getop('matrixelement', z.op, id='test_comparisons'))
        s = SL(2, z, getop('matrixelement', z.op, id='test_comparisons'))
        g = GL(2, z, getop('matrixelement', z.op, id='test_comparisons'))
        self.assertEqual(g<m, True)
        self.assertEqual(s<g, True)
        self.assertEqual(s<m, True)
        self.assertEqual(s<=g, True)
        self.assertEqual(s>g, False)
        self.assertEqual(g<s, False)
        self.assertEqual(g==s, False)
        self.assertEqual(m==g, False)
        self.assertEqual(m>s, True)
        self.assertEqual(g>s, True)
        self.assertEqual(g>=s, True)

    def test_permutations(self):
        self.assertEqual(Permutation([1, 3, 4])[3], 4)
        self.assertEqual(Permutation([1, 3, 4])[4], 1)
        self.assertEqual(Permutation([1, 3]),
            Permutation([3, 1]))
        self.assertEqual(Permutation([1, 3]),
            Permutation([2, 1])*
            Permutation([2, 3])*
            Permutation([2, 1]))
        self.assertEqual(Permutation([1, 3]).sign(),-1)
        self.assertEqual(Permutation([1, 3, 5]).sign(),1)
        self.assertEqual((Permutation([3, 5])*
            Permutation([6, 4, 7])*
            Permutation([1, 12, 14])*
            Permutation([5, 3])).sign(),1)
        self.assertEqual((Permutation([3, 5])*
            Permutation([5, 3])).sign(),1)

    def test_permutationgroup(self):
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

    def test_generatorgroup(self):
        for i in range(5,20):
            self.assertEqual(Z(i),
                GeneratorGroup(1 if i>1 else (0 if i==1 else None), getop('addmod', i)))
        self.assertEqual(PermutationGeneratorGroup([[[1,2,3]],[[1,2]]]), S(3))

    def test_permutationmatrix(self):
        self.assertEqual(Aff(D(2)), GeneratorGroup(Matrix([
            [Permutation([1,2]),Permutation([1,2])],
            [0,1]]), Aff(D(2)).op))
        self.assertEqual(Aff(D(2)), GeneratorGroup(Matrix([
            [Permutation([1,2]),Permutation([])],
            [0,1]]), Aff(D(2)).op))

    def test_quaternion(self):
        self.assertEqual(Quaternion(0,1,0,0)*Quaternion(0,1,0,0),-1)
        self.assertEqual(Quaternion(0,0,1,0)*Quaternion(0,0,1,0),-1)
        self.assertEqual(Quaternion(0,0,0,1)*Quaternion(0,0,0,1),-1)
        self.assertEqual(Quaternion(0,1,0,0)*Quaternion(0,0,1,0)*Quaternion(0,0,0,1),-1)
        self.assertEqual(Q(8), Dic(2))
        for i in range(1,20):
            self.assertEqual(len(Dic(i)), 4*i)

if __name__ == '__main__':
    unittest.main()
