from sage.rings.integer import GCD_list
from sage.geometry.polyhedron.constructor import Polyhedron
import itertools
from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
import numpy as np
from pprint import pprint
#from heapq import heapify, heappush, heappop

from sage.quadratic_forms.qfsolve import qfsolve, qfparam
from sympy.solvers.diophantine import *
import sympy


#local imports:
import coxiter

# taken from https://gist.github.com/sirodoht/ee2abe82eca70f5b1869
from operator import mul, mod

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)


def crt(m, a):
    M = reduce(mul, m) # the product of m elements
    m_i = [M / item for item in m]
    b = map(mod, m_i, m)
    g, k, l = map(egcd, b, m)
    g, k, l = zip(g, k, l) # transpose g, k and l arrays
    t = map(mod, k, m)
    e = map(mul, m_i, t)
    
    x_sum = sum(map(mul, a, e))
    x = x_sum % M
    return x



def qform(B):
    C = 2*B
    Q = QuadraticForm(QQ, C)
    return Q


def solve_mod_primes(n,diff,primes):
    def solve_mod_p(n,diff,p):
        Box = itertools.product(*[range(p) for i in range(n)])
        return [x for x in Box if diff(x)%p == 0]
    
    def vectors_multiply(vectors, primes):
        n = len(vectors[0])
        return tuple(crt(primes,[v[i] for v in vectors]) for i in range(n))
    sols_mod_primes = [solve_mod_p(n,diff,p) for p in primes]
    sol_combinations = [i for i in itertools.product(*sols_mod_primes)]
    merged_solutions = set([ vectors_multiply(sol, primes) for sol in sol_combinations]) # need to filter repeated vectors
    modulus = prod(primes) # m=modulus; plan: rebalance to [-m/2..m/2], iterate y over Box/m and look for Q(my+x)=k. Filter y's somehow.
    return merged_solutions, modulus

def define_primes(boundary):
    return [2,3,5]

def Solve_equation(n, boundary, diff):
    remainders, modulus = solve_mod_primes(n,diff, define_primes(boundary))
    boundary = round(boundary/modulus) + 1
    BOX = itertools.product(*[range(-boundary,boundary) for i in range(n)])
    sols = []
    for y in BOX:
        for x in remainders:
            v = vector(y)*modulus+vector(x)
            if (diff(v) == 0):
                sols.append(v)
    return sols



def ParallelepipedContains(P,v):#v is a row, P is a matrix with polytope generators as rows
    Q = matrix(v)*P.inverse()
    return all( (c < 1) and (c>=0) for c in Q.list())

def BoundingBox(vectors):
    n = len(vectors[0])
    negative = [sum(v[i] for v in vectors if v[i]<0) for i in range(n)]
    positive = [sum(v[i] for v in vectors if v[i]>0) for i in range(n)]
    return itertools.product(*[range(negative[i],positive[i]+1) for i in range(len(negative))])

def GetIntegerPoints(m):
    return [v for v in BoundingBox(m.rows()) if ParallelepipedContains(m,v)]


def ReduceToV1(V1, u, k): # reduces Q(v1+u)=k to Q(v1+u1)=k1, where u1 is in V1
    E = V1.vector_space()
    u1 = sum(e.inner_product(u)*eps for (e, eps) in itertools.collate(E.basis(), E.dual_basis()))
    k1 = k - u.inner_product(u) + u1.inner_product(u1)
    return u1, k1






def NegativeVector(V):
    D, T = QuadraticForm(QQ, V.inner_product_matrix()).rational_diagonal_form(True)
    D = D.matrix()
    #print('\n'.join(["M=",repr(M),"D=",repr(D),"T=",repr(T)]))
    #print [(D[i][i], T.column(i)) for i in range(D.nrows())]
    vectors = [T.column(i) for i in range(D.nrows()) if D[i][i]<0]
    assert len(vectors) == 1
    v0 = vectors[0]
    v0 = v0 / gcd(v0.list())
    return V(v0)





            
class VinAl:
    def __init__(s, M, v0=None):
        s.M = M
        s.n = s.M.ncols()  # n-1 = dimension of hyperbolic space
        s.V = FreeModule(ZZ,s.n, inner_product_matrix=s.M) # created a quadratic lattice
        if v0 is None:
            s.v0 = NegativeVector(s.V)
        else:
            assert s.n == len(v0)
            s.v0 = s.V(v0)
            
        s.V1 = s.V.submodule(matrix([s.v0.dot_product(m) for m in s.M.columns()]).right_kernel()) # V1 = <v0>^\perp
        s.M1 = s.V1.gram_matrix()
        s.mu = (sorted(s.M1.eigenvalues()))[0]
        
        s.W=[s.V(w) for w in GetIntegerPoints(s.V1.matrix().insert_row(s.n-1, s.v0))] # w_1..w_m
        s.W.sort(key = lambda x: -s.v0.inner_product(x))

        s.En=abs(s.M.det()/GCD_list(s.M.adjoint().list()))
        s.root_lengths = [k for k in range(1,2*s.En+1) if ((2*s.En)%k == 0)]

        #print([v0.inner_product(w) for w in W])
        s.check_validity()
        s.roots = []
        print("Vinberg algorithm initialized\n")

    def Print(s): # change to __str__ and/or __repr__
        print("W:")
        print(s.W)
        print("V1:")
        print(s.V1.vector_space())
        print(s.V1.gram_matrix())
        
    def check_validity(s):
        assert s.M.is_square()
        assert s.M.is_symmetric()
        assert s.v0.inner_product(s.v0) < 0 # checking <v0,v0> < 0
        assert s.V1.gram_matrix().is_positive_definite()
        assert all(0>=s.v0.inner_product(w) and s.v0.inner_product(w)>s.v0.inner_product(s.v0) for w in s.W)
        #assert len(s.W) == s.V.submodule(s.V1.gens()+(v0,)).index_in(s.V) # strange error

    def IsRoot(s, v):
        return all( (2*v.inner_product(e))%v.inner_product(v)==0 for e in s.V.gens() )

    def IsNewRoot(s, v):
        if s.V.are_linearly_dependent(s.roots[:s.V.degree()-1]+[v,]):
            return False
        return s.IsRoot(v) and all(v.inner_product(root)<=0 for root in s.roots)


    def is_FundPoly(s):
            if len(s.roots)<1:
                return False
            M = [[ t.inner_product(r) for t in s.roots] for r in s.roots]
            print(M)
            return coxiter.run(M, s.n)
        
    def FindRoots(s):
        s.roots = s.FundCone()
        for root in s.NextRoot():
            s.roots.append(root)
            if s.is_FundPoly():
                print('Fundamental Polyhedron constructed, roots:')
                print(s.roots)
                return

    def FundCone(s):
        V1_roots = [v for k in s.root_lengths for v in s.Roots_decomposed_into(s.V([0]*s.n), k) if s.IsRoot(v)]
        print('roots in V1:', V1_roots)
        cone = Cone([[0]*s.n]).dual()
        #print('cone', cone.rays())
        for root in V1_roots:
            halfplane = Cone([root]).dual()
            #print('halfplane', halfplane.rays())
            if cone.intersection(halfplane).dim() == s.n:
                cone = cone.intersection(halfplane)
            else:
                cone = cone.intersection(Cone([-root]).dual())
        #print('cone', cone.rays())
        print('FundCone returned')
        return [s.V(r) for r in cone.dual().rays()]
        
        
    def Roots_decomposed_into(s, a, k): #k is desired inner square, a is a non-V1 component
        boundary = round(sqrt(k/s.mu))+max(abs(a[i]) for i in range(s.n))
        def V1_vector(x):
            return sum(x[i]*s.V1.gens()[i] for i in range(s.n-1))
        def diff(x):
            v1 = V1_vector(x)
            return (a+v1).inner_product(a+v1)-k
        return [V1_vector(x)+a for x in Solve_equation(s.n-1, boundary, diff)]

    def IterateRootDecompositions(s, stop=-1): # iterates pairs (w_i + c v_0, ||a||) from minimum, infinity or `stop` times
        candidates = {k:1 for k in s.root_lengths} # dictionary of vector numbers for each k: we have a series of vectors {W}\0, {v0 + W}, etc. The number is the position in this series.
        def cand_a(n): # the non-V1 component of vector under number n
            return s.W[n%len(s.W)]+ (n//len(s.W))*s.v0
        while True:
            k = min(s.root_lengths, key=lambda k: -s.v0.inner_product(cand_a(candidates[k]))/math.sqrt(k)) # can be optimized with heapq
            #print 'IterateDecompositions returns ', cand_a(candidates[k]), k
            yield cand_a(candidates[k]), k
            candidates[k]+=1
            stop-=1
            if stop == 0:
                return

    def NextRoot(s):
      for a, k in s.IterateRootDecompositions():
        print(a, k, -a.inner_product(s.v0)/math.sqrt(k))
        new_roots = [v for v in s.Roots_decomposed_into(a, k) if s.IsNewRoot(v)]
        if len(new_roots)>0:
            print 'root candidates', new_roots
        for root in new_roots:
            yield root
            

# M is an inner product (quadratic form), v0 is a chosen vector
M = diagonal_matrix(ZZ,[-1,1,1])
v0 = [1,0,0,0]
A = VinAl(M)
