import timeit
#from qsolve import squares_sum_solve as solve
from vinal import VinAl
#from inf_symm import *
from forms import *
import tests
import unittest

import PyNormaliz
C = PyNormaliz.Cone(cone = [[1,0],[0,1]])


def test_vinal():
  A = VinAl([[0,1,0],[1,0,0],[0,0,1]])
  print(A.n)
  print(rational_diagonal_form(A.Q))
  print(A.En)


def test_forms():
  print(parallelepiped_integer_points(Matrix([[2,1],[1,2]])))

def test_symm():
  print (inner(roots[0],roots[1]))



test_vinal()
#test_symm()
#test_forms()


suite = unittest.TestLoader().loadTestsFromModule(tests)
unittest.TextTestRunner(verbosity=2).run(suite)
