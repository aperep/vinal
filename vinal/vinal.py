from sympy import *
from cached_property import cached_property
import math
import functools
from fractions import Fraction

from forms import rational_diagonal_form
from cone import Cone
from qsolve import squares_sum_solve
import coxiter

class Lattice: 
  def __init__(self, Q, v0=None):
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
    denominators = [x.q for x in tuple(vector)]
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
 

class VinAl(Lattice):
  '''
  We always work in diagonal basis, except input and output
  '''
  def __init__(self, Q, v0=None):
    super().__init__(Q,v0)
    self.roots = None
    self.roots_diag = None

  @cached_property
  def root_lengths(self): # possible lengths of roots
    return [k for k in range(1,2*self.En+1) if ((2*self.En)%k == 0)]

  def is_root(self, v, diag = False):
    v_length = self.dot(v,v, diag=diag)
    assert v_length != 0
    return all( 2*d % v_length == 0 for d in self.dot(v, self.basis(diag=diag), diag=diag) ) 

  def is_new_root(self, v, diag = False):
    roots = self.roots_diag if diag else self.roots
    M = roots[:, :self.n-1].col_insert(0, v)
    if M.rank() < M.cols:
      return False
    return self.is_root(v, diag=diag) and all( dot <= 0 for dot in self.dot(v,roots, diag=diag) )

  def root_types(self, diag = False, stop=-1): # iterates pairs (a = w_i + c v_0, ||a||) from minimum, infinity or `stop` times
    a_num = {k:1 for k in self.root_lengths} # each possible length k we store an index of a from series w_1, ..., v0, v0+w_1, ... 
    W = self.W(diag=diag)
    v0 = Matrix([1]+[0]*(n-1)) if diag else self.v0
    def a(k): # the non-V1 component of vector a for length k
      n = a_num[k]
      m = len(W)
      return W[n%m] + (n//m)*v0
    while True:
      k = min(self.root_lengths, key=lambda k: -self.dot(v0,a(k),diag=diag)/math.sqrt(k)) 
      yield a(k), k
      a_num[k]+=1
      stop-=1
      if stop == 0:
        return

  def add_root(self, root, diag = False):
    root_nondiag = self.from_diag(root) if diag else root
    root_diag = root if diag else self.to_diag(root)
    self.roots = self.roots.col_insert(self.roots.cols, root_nondiag)
    self.roots_diag = self.roots.col_insert(self.roots.cols, root_diag)


  def roots_of_type(self, a, k, diag = False): #k is desired length squared, a is a non-V1 component
    return [v if diag else self.from_diag(v) for s in self.shifts_in_V1 for v in self.roots_with_shift(a+s, k, diag=True)]

  def roots_with_shift(self, a, k, diag = False): # TODO: we should first find roots for a given shift and then take all V1-shifts
    '''
        Here we solve the equation (a+v1, a+v1) == k for a vector v1 in V1.
    ''' 
    if not diag:   
      a = self.to_diag(a)
    q = [self.Q_diag[i,i] for i in range(1, self.n)]
    c = k - self.Q_diag[0,0]*a[0,0]**2
    for solution in squares_sum_solve(q, c, offset = list(a)[1:]):
      result = Matrix([0]+solution) + Matrix(a)
      yield result if diag else self.from_diag(result) # TODO: check correctness

  @functools.lru_cache()
  def roots_in_v0_perp(self, diag = False):
    return [v for k in self.root_lengths 
      for v in self.roots_of_type(zeros(self.n,1), k, diag=diag) 
        if self.is_root(v, diag=diag)]

  def next_root(self, diag = False,  stop = -1):
    for a, k in self.root_types(diag=diag, stop=stop):
      print('next_root: working with type ', a, k)#, -a.inner_product(s.v0)/math.sqrt(k))
      new_roots = [v for v in self.roots_of_type(a, k, diag=diag) if self.is_new_root(v, diag=diag)]
      if len(new_roots)>0:
        print(f'new root candidates {"in diagonal basis" if diag else ""}', new_roots)
        pass
      for root in new_roots:
        print('next_root: yielding root ', root)
        yield root # what exactly we want to do here?


  def run(self, stop = -1):
    self.roots = self.fundamental_cone() # roots are in diagonal coordinates !! 
    self.roots_diag = self.to_diag(self.roots)
    if not self.finished():
      for root in self.next_root(stop=stop):
        self.add_root(root, diag=False)
        print('roots found: {0}, they are:\n{1}'.format(self.roots.cols,self.roots))
        if self.finished():
          break
    print('Fundamental Polyhedron constructed, roots:')
    print(self.roots) # TODO: return to non-diagonal coordinates
    return self.roots.T 



  def fundamental_cone(self, diag = False): 
    cone = Cone([[0]*self.n])
    for root in self.roots_in_v0_perp(diag=False):
      if cone.intersects(root):
        cone.append(root)
    print('FundCone constructed, roots: ',cone.rays)
    rays = Matrix(list(cone.rays)).T
    return self.to_diag(rays) if diag else rays


  def finished(self, diag = False):
    roots = self.roots_diag if diag else self.roots
    if len(roots)<2:
      return False
    M = self.dot(roots, roots, diag=diag)
    print('checking polyhedron with Gram matrix:\n',M)
    return coxiter.run(M.tolist(), self.n)


if __name__ == '__main__':
  A = [[0,1,0],[1,0,0],[0,0,1]]
  V = VinAl(A)
  roots = V.run()
  print(roots)
  D = diag(-3,5,1,1)
  VD = VinAl(D)
