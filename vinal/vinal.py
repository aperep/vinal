#!/usr/bin/env sage -python

from sympy import *
from cached_property import cached_property
import math
import functools
from fractions import Fraction

from lattice import Lattice
from cone import Cone
from qsolve import squares_sum_solve
import coxiter
 
class VinAl(Lattice):
  '''
  We always work in diagonal basis, except input and output
  '''
  def __init__(self, Q, v0=None, blocks = 'sympy'):
    super().__init__(Q,v0)
    self.roots = None
    self.roots_diag = None
    if blocks == 'sage':
      import blocks_sage_original as blocks
      # substituting methods as in https://stackoverflow.com/a/2982/7626757 (may use types.MethodType as well)
      self.roots_of_type = blocks.roots_of_type.__get__(self) 
      self.root_types = blocks.root_types.__get__(self)
      self.fundamental_cone = blocks.fundamental_cone.__get__(self)

  @cached_property
  def root_lengths(self): # possible lengths of roots
    return [k for k in range(1,2*self.En+1) if ((2*self.En)%k == 0)]

  def is_root(self, v, diag = False):
    v_length = self.dot(v,v, diag=diag)
    assert v_length != 0
    return all( 2*d % v_length == 0 for d in self.dot(v, self.basis(diag=diag), diag=diag) ) 

  def is_new_root(self, v, diag = False):
    if not self.is_root(v, diag=diag):
      return False
    roots = self.roots_diag if diag else self.roots
    M = roots[:, :self.n-1].col_insert(0, v)
    if M.rank() < M.cols:
      return False
    return all( dot <= 0 for dot in self.dot_matrix(v,roots, diag=diag) )

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
      print(f'next_root: working with type a={a}, k={k} ', a, k)#, -a.inner_product(s.v0)/math.sqrt(k))
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
    if roots.cols<2:
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
