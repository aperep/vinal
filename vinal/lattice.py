from sympy import *
from cached_property import cached_property
import math
import functools
from fractions import Fraction

from forms import rational_diagonal_form

class Lattice: 
  def __init__(self, Q, v0=None, d=None):
    '''
    here we establish necessary coordinate systems
    Q - integer-valued quadratic form on Z^n in a form of an array of arrays of integers
    v0 - optional starting negative vector
    '''
    self.Q = Matrix(Q) # quadratic form on the ambient space
    self.n = self.Q.rows # dimension of ambient space

    if v0 == None:
      self.Q_diag, self.basis_diag  = rational_diagonal_form(self.Q)
      self.v0                       = self.basis_diag[:,0] # do we need that?
      self.V1_basis_diag            = self.basis_diag[:,1:]
      self.basis_diag_inverse       = self.basis_diag ** -1
    else:
      raise NotImplementedError # TODO: First cover with tests, then fix this case
      self.v0                       = Matrix(v0)
      V1_skew_basis                 = (self.Q*self.v0).T.nullspace()
      Q1                            = V1_skew_basis.T * self.Q * V1_skew_basis
      Q1_diag, self.V1_basis        = rational_diagonal_form(Q1)
      self.basis_diag               = [self.v0] + self.V1_basis # TODO: inconsistent dimension
      self.Q_diag                   = diag(self.v0.T * self.Q * self.v0, Q1_diag)
    
    assert self.dot(self.v0, self.v0) < 0

  def basis(self, diag = False):
    return self.basis_diag_inverse if diag else eye(self.n) 

  @cached_property
  def En(self): # higher ???, which is an upper bound of root length 
    adjoint_elements = [self.Q.cofactor(i,j) for i in range(self.n) for j in range(self.n)]
    return abs(self.Q.det()/gcd(adjoint_elements))

  def to_diag(self, v):     return self.basis_diag_inverse * v
  def from_diag(self, v):   return self.basis_diag * v

  def dot(self, a, b, diag = False):
    Q = self.Q_diag if diag else self.Q
    product = (a.T*Q*b)
    if product.shape == (1, 1):
      return product[0,0]
    else:
      return product
  
  def vectors(self, M):
    return [M.col(q) for q in range(M.cols)]

  def to_unit_cube(self, vector): # returns a shift of vector in diagonal coords by diag_lattice in the unit cube
    return vector.applyfunc(lambda x: x%1)
  
  def order(self, vector):
    '''
    vector is an element of lattice (in diagonal coordinates?)
    Returns order of vector image in the quotient lattice
    '''
    denominators = [int(x) for x in tuple(vector)]
    return functools.reduce(lambda x,y:x*y//math.gcd(x,y),denominators)

  @cached_property 
  def shifts(self): # all shifts (in unit cube) of diagonal lattice that comprise the whole lattice
    multiples = lambda vector: {ImmutableMatrix(i*vector) for i in range(self.order(vector))}
    basis_multiples = [multiples(v) for v in self.vectors(self.basis(diag=True))]
    minkowski_sum = lambda s1,s2: {ImmutableMatrix(self.to_unit_cube(v1+v2)) for v1 in s1 for v2 in s2}
    return functools.reduce( minkowski_sum, basis_multiples)

  @cached_property 
  def shifts_in_V1(self): # all shifts that lie in V1
    return {s for s in self.shifts if s[0,0]==0}

  @cached_property 
  def W_diag(self): # all shifts of V1+<v0> that give Z^n
    first_coords = [s[0,0] for s in self.shifts if s[0,0]!=0]
    if len(first_coords)==0:
      return [zeros(self.n,1)]
    primitive_vector = [s for s in self.shifts if s[0,0]==min(first_coords)][0]
    return [self.to_unit_cube(i*primitive_vector) for i in range(self.order(primitive_vector))]
 
  def W(self, diag = False): # all shifts of V1+<v0> that give Z^n
    return self.W_diag if diag else [self.from_diag(w) for w in self.W_diag]
 
