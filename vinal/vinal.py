from sympy import *
from cached_property import cached_property

from forms import *

class Lattice: 
  def __init__(self, Q, v0=None):
    '''
    here we establish necessary coordinate systems
    Q - integer-valued quadratic form on Z^n in a form of an array of arrays of integers
    keyword parameters:
    v0 - optional starting negative vector
    '''
    self.Q = Matrix(Q) # quadratic form on the ambient space
    self.n = self.Q.shape[0] # dimension of ambient space
    self.basis = [eye(self.n).col(i) for i in range(self.n)]

    if v0 == None:
      self.Q_diag, self.basis_diag = rational_diagonal_form(self.Q)
      self.v0 = self.basis_diag[:,0] # do we need that?
      self.V1_basis = self.basis_diag[:,1:]
    else:
      self.v0 = Matrix(v0)
      V1_skew_basis = (self.Q*self.v0).T.nullspace()
      Q1 = V1_skew_basis.T * self.Q * V1_skew_basis
      Q1_diag, self.V1_basis = rational_diagonal_form(Q1)
      self.basis_diag = [self.v0] + self.V1_basis
      self.Q_diag = diag(self.v0.T * self.Q * self.v0, Q1_diag)
      
    assert (self.v0.T * self.Q * self.v0)[0,0] < 0

  @cached_property 
  def W(self): # all shifts of orthogonal lattice that give Z^n
    return sorted(parallelepiped_integer_points(self.basis_diag), key = lambda x: -(self.v0.T * self.Q * x)[0,0]) # 
 
  @cached_property
  def En(self): # higher ???, which is an upper bound of root length 
    adjoint_elements = [self.Q.cofactor(i,j) for i in range(self.n) for j in range(self.n)]
    return abs(self.Q.det()/gcd(adjoint_elements))




class VinAl(Lattice):
  def __init__(self, Q, v0=None):
    super().__init__(Q,v0)
    self.roots = []


  @cached_property
  def root_lengths(self): # possible lengths of roots
    return [k for k in range(1,2*self.En+1) if ((2*self.En)%k == 0)]

  def is_root(self, v):
        return all( ( (2*v.T*self.Q*e) % (v.T*self.Q*v) ) == 0 for e in self.basis )

  def is_new_root(self, v):
        vector_system = Matrix(self.roots[:self.n-1]+[v,]).T
        if vector_system.rank()< vector_system.cols:
            return False
        return self.is_root(v) and all(v.T*self.Q*root <=0 for root in self.roots)



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
          yield Matrix([0]+solution) + a


  @cached_property
  def roots_in_v0_perp(self): # possible lengths of roots
      return [v for k in self.root_lengths for v in self.roots_of_type([0]*self.n, k) if self.is_root(v)]

  def next_root(self):
      for a, k in self.root_types():
          #print(a, k, -a.inner_product(s.v0)/math.sqrt(k))
        new_roots = [v for v in self.roots_of_type(a, k) if self.is_new_root(v)]
        if len(new_roots)>0:
            #print('new root candidates', new_roots)
            pass
        for root in new_roots:
            yield root # what exactly we want to do here?


  def run(self):
        self.roots = self.fundamental_cone()
        for root in self.next_root():
            self.roots.append(root)
            print('roots found: {0}, they are:\n{1}'.format(len(self.roots),self.roots))
            if self.finished():
                print('Fundamental Polyhedron constructed, roots:')
                print(self.roots)
                return


  def fundamental_cone(self): # we don't have polyhedral cones in sympy, need to rethink.
  # either take algorithm from sagemath, or use something like PyNormaliz, pyparma, or pycddlib - they implement double description of polyhedra;
  # for self-contained implementation, look at page 5-8 of http://www.lab2.kuis.kyoto-u.ac.jp/~avis/courses/pc/2010/notes/lec5.pdf or function `initial_pair` at https://github.com/sagemath/sage/blob/e8633b09919542a65e7e990c8369fee30c7edefd/src/sage/geometry/polyhedron/double_description.py
  # Trying to implement this version ftp://ftp.math.tu-berlin.de/pub/Lehre/LinOpt/WS06/Uebung15.11.06/fukuda96double.pdf
       
        for root in self.roots_in_v0_perp:
            halfplane = Cone([root]).dual()
            #print('halfplane', halfplane.rays())
            if cone.intersection(halfplane).dim() == s.n:
                cone = cone.intersection(halfplane)
            else:
                cone = cone.intersection(Cone([-root]).dual())
        #print('cone', cone.rays())
        print('FundCone returned',[s.V(r) for r in cone.dual().rays()])
        return [s.V(r) for r in cone.dual().rays()]


  def finished(self):
            if len(s.roots)<1:
                return False
            M = [[ t.inner_product(r) for t in s.roots] for r in self.roots]
            print('checking polyhedron with Gram matrix')
            print(Matrix(M))
            return coxiter.run(M, self.n)
