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

print (3/2)%5

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

def Get_Iterate_ak():
    pass
#

def IterateDecompositions(kroot_lengths,  W, v0, stop=-1): # iterates pairs (w_i + c v_0, ||a||) from minimum, infinity or `stop` times
    #print kroot_lengths
        candidates = {k:0 for k in kroot_lengths} # dictionary of vector numbers for each k: we have a series of vectors {W}, {v0 + W}, etc. The number is the position in this series.
        def cand_a(n): # the non-V1 component of vector under number n
            return W[n%len(W)]+ (n//len(W))*v0

        while True:
            k = min(kroot_lengths, key=lambda k: -v0.inner_product(cand_a(candidates[k]))/math.sqrt(k)) # can be optimized with heapq
            #print 'IterateDecompositions returns ', cand_a(candidates[k]), k
            yield cand_a(candidates[k]), k
            candidates[k]+=1
            stop-=1
            if stop == 0:
                return

def ReduceToV1(V1, u, k): # reduces Q(v1+u)=k to Q(v1+u1)=k1, where u1 is in V1
    E = V1.vector_space()
    u1 = sum(e.inner_product(u)*eps for (e, eps) in itertools.collate(E.basis(), E.dual_basis()))
    k1 = k - u.inner_product(u) + u1.inner_product(u1)
    return u1, k1

def v1Box(u1, k1): # get bounding box for v1
    pass

def Roots(V, V1, a, k): #k is desired inner square, a is a non-V1 component
    M1 = V1.gram_matrix()
    #print M1
    n = V1.degree()
    #print('n',n)
    mu = (sorted(M1.eigenvalues()))[0]
    boundary = round(sqrt(k/mu))+max(abs(a[i]) for i in range(n))
    def V1_vector(x):
        return sum(x[i]*V1.gens()[i] for i in range(n-1))
    def diff(x):
        v1 = V1_vector(x)
        #print('v1', x, v1)
        #print(v1, a)
        return (a+v1).inner_product(a+v1)-k
    return [V1_vector(x)+a for x in Solve_equation(n-1, boundary, diff)]
#return (V1(v) for v in BoundingBox([a*k for a in V1.gens()]+[-a*k for a in V1.gens()]) if check(v))
# google quadraticform .short_primitive_vector_list_up_to_length, .solve(k), .vectors_by_length
# modulo p reduction!

def GetRoot(V, V1, M,W,v0,  roots):
    def IsRoot(v):
        k = v.inner_product(v)
        if V.are_linearly_dependent(roots+[v,]) and (len(roots) < V.degree()):
            return False
        return all( (2*v.inner_product(e))%k==0 for e in V.gens() ) and all(v.inner_product(root)<=0 for root in roots)
    En=abs(M.det()/GCD_list(M.adjoint().list()))
    kroot_lengths = [k for k in range(1,2*En+1) if ((2*En)%k == 0)]
    for a, k in IterateDecompositions(kroot_lengths, W, v0):
        print(a, k, -a.inner_product(v0)/math.sqrt(k))
        #print 'solutions', Roots(V, V1, a, k)
        new_roots = [v for v in Roots(V, V1, a, k) if IsRoot(v)]
        if len(new_roots)>0:
            print 'root candidates', new_roots
            return new_roots[0]



def vinberg_algorithm(M, v0=None): # M is an inner product (quadratic form), v0 is a chosen vector
        #pprint("Vinberg algorithm started\n")
        print("Vinberg algorithm started\n")
        assert M.is_square()
        assert M.is_symmetric()
        n = M.ncols()  # n-1 = dimension of hyperbolic space
        V = FreeModule(ZZ,n, inner_product_matrix=M) # created a quadratic lattice
        if v0 is None:
            v0 = NegativeVector(V)
        else:
            assert n == len(v0)
            v0 = V(v0)
        assert v0.inner_product(v0) < 0 # checking <v0,v0> < 0
        V1 = V.submodule(matrix([v0.dot_product(m) for m in M.columns()]).right_kernel()) # V1 = <v0>^\perp
        assert V1.gram_matrix().is_positive_definite()
        print(V1.gram_matrix())
        W=[V(w) for w in GetIntegerPoints(V1.matrix().insert_row(n-1, v0))] # w_1..w_m
        W.sort(key = lambda x: -v0.inner_product(x))
        print("W:")
        print(W)
        #print([v0.inner_product(w) for w in W])
        Vprime=V.submodule(V1.gens()+(v0,)) # V'= V1 + <v0>
        assert len(W) == Vprime.index_in(V)
        assert all(0>=v0.inner_product(w) and v0.inner_product(w)>v0.inner_product(v0) for w in W)
        print("V1:")
        print(V1.vector_space())
        roots = []
        while not is_FundPoly(roots):
            roots.append(GetRoot(V,V1,M,W,v0, roots))
            print 'Roots:', roots

def is_FundPoly(roots):
    return len(roots)>3

vinberg_algorithm(diagonal_matrix(ZZ,[-3,5,1,1]), [1,0,0,0])



