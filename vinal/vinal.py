from sympy import *
from cached_property import cached_property
import math
import functools
from fractions import Fraction

from forms import rational_diagonal_form, parallelepiped_integer_points
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
    self.basis = eye(self.n)

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

  @cached_property 
  def W(self): # all shifts of orthogonal lattice that give Z^n
    W = [self.primitive_vector]
    while True:
      vector = self.to_unit_cube(W[-1]+self.primitive_vector)
      if vector[0,0]==0:
        break
      W.append(vector)
    W = [self.from_diag(w) for w in W]
    return sorted(W, key = lambda x: -self.dot(self.v0, x)) # 
 
  @cached_property
  def En(self): # higher ???, which is an upper bound of root length 
    adjoint_elements = [self.Q.cofactor(i,j) for i in range(self.n) for j in range(self.n)]
    return abs(self.Q.det()/gcd(adjoint_elements))

  def dot_Q(self, a, b, Q):
    product = (a.T*self.Q*b)
    if product.shape == (1, 1):
      return product[0,0]
    else:
      return product
  
  def dot(self, a, b):
    return self.dot_Q(a,b,self.Q)

  def dot_diag(self, a, b):
    return self.dot_Q(a,b,self.Q_diag)

  def to_diag(self, v):
    return self.basis_diag_inverse * v

  def from_diag(self, v):
    return self.basis_diag * v

  def vectors(self, M):
    return [M.col(q) for q in range(M.cols)]

  def to_unit_cube(self, vector): # returns a shift of vector in diagonal coords by diag_lattice in the unit cube
    return vector.applyfunc(lambda x: x%1)
  
  def order(self, vector):
    denominators = [x.denominator for x in tuple(vector) if type(x)==Fraction]
    if len(denominators)==0:
      return 1
    return functools.reduce(lambda x,y:x*y//math.gcd(x,y),denominators)

  @cached_property 
  def primitive_vector(self):    
    def common_generator(v1,v2):
      order1 = self.order(v1)
      order2 = self.order(v2)
      if order1 % order2 == 0:
        return self.to_unit_cube(v1)
      if order2 % order1 == 0:
        return self.to_unit_cube(v2)
      return self.to_unit_cube(v1+v2)
    return  functools.reduce(common_generator,self.vectors(self.basis_diag_inverse))

class VinAl(Lattice):
  def __init__(self, Q, v0=None):
    super().__init__(Q,v0)
    self.roots = None

  @cached_property
  def root_lengths(self): # possible lengths of roots
    return [k for k in range(1,2*self.En+1) if ((2*self.En)%k == 0)]

  def is_root(self, v):
    v_length = self.dot(v,v)
    assert v_length != 0
    return all( 2*d % v_length == 0 for d in self.dot(v, self.basis) ) 

  def is_new_root(self, v):
        M = self.roots[:, :self.n-1].col_insert(0, v)
        if M.rank() < M.cols:
            return False
        return self.is_root(v) and all( dot <= 0 for dot in self.dot(v,self.roots) )



  def root_types(self, stop=-1): # iterates pairs (a = w_i + c v_0, ||a||) from minimum, infinity or `stop` times
        a_num = {k:1 for k in self.root_lengths} # each possible length k we store an index of a from series w_1, ..., v0, v0+w_1, ... 
        def a(k): # the non-V1 component of vector a for length k
            n = a_num[k]
            m = len(self.W)
            return self.W[n%m] + (n//m)*self.v0
        while True:
            k = min(self.root_lengths, key=lambda k: -(self.v0.T*self.Q*a(k))[0,0]/math.sqrt(k)) 
            yield a(k), k
            a_num[k]+=1
            stop-=1
            if stop == 0:
                return


  def roots_of_type(self, a, k): #k is desired length squared, a is a non-V1 component
        '''
        Here we solve the equation (a+v1, a+v1) == k for a vector v1 in V1.
        '''
        q = [self.Q_diag[i,i] for i in range(1, self.n)]
        c = k - self.Q_diag[0,0]*a[0]**2
        for solution in squares_sum_solve(q, c, offset = a[1:]):
          yield self.from_diag( Matrix([0]+solution) + Matrix(a) ) # TODO: check correctness


  @cached_property
  def roots_in_v0_perp(self): # possible lengths of roots
    #print([v for v in self.roots_of_type([0]*self.n, 2)])
    return [v for k in self.root_lengths for v in self.roots_of_type([0]*self.n, k) if self.is_root(v)]

  def next_root(self):
      for a, k in self.root_types():
        print('next_root: working with type ', a, k)#, -a.inner_product(s.v0)/math.sqrt(k))
        new_roots = [v for v in self.roots_of_type(a, k) if self.is_new_root(v)]
        if len(new_roots)>0:
            print('new root candidates', new_roots)
            pass
        for root in new_roots:
            print('next_root: yielding root ', root)
            yield root # what exactly we want to do here?


  def run(self):
        self.roots = self.fundamental_cone() # roots are in diagonal coordinates !! 
        if not self.finished():
          for root in self.next_root():
            self.roots = self.roots.col_insert(self.roots.cols, root)
            print('roots found: {0}, they are:\n{1}'.format(self.roots.cols,self.roots))
            if self.finished():
                break
        print('Fundamental Polyhedron constructed, roots:')
        print(self.roots) # TODO: return to non-diagonal coordinates
        return self.roots.T 



  def fundamental_cone(self): 
    cone = Cone([[0]*self.n])
    for root in self.roots_in_v0_perp:
        if cone.intersects(root):
            cone.append(root)
    print('FundCone constructed, roots:',cone.rays)
    return Matrix(list(cone.rays)).T


  def finished(self):
            if len(self.roots)<2:
                return False
            M = self.dot(self.roots, self.roots)
            print('checking polyhedron with Gram matrix')
            print(M)
            return coxiter.run(M.tolist(), self.n)


if __name__ == '__main__':
  A = [[0,1,0],[1,0,0],[0,0,1]]
  V = VinAl(A)
  roots = V.run()
  print(roots)
  D = diag(-3,5,1,1)
  VD = VinAl(D)
