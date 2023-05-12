from sage.all_cmdline import *   # import sage library

from sage.rings.integer import GCD_list
import itertools
import numpy as np

import qsolve
import coxiter



def GetIntegerPoints(m):
    n = len(m.rows()[0])
    negative = [sum(v[i] for v in m.rows() if v[i]<0) for i in range(n)]
    positive = [sum(v[i] for v in m.rows() if v[i]>0) for i in range(n)]
    BoundingBox = itertools.product(*[range(negative[i],positive[i]+1) for i in range(len(negative))])
    def ParallelepipedContains(v):#v is a row, m is a matrix with polytope generators as rows
        Q = matrix(v)*m.inverse()
        return all( (c < 1) and (c>=0) for c in Q.list())
    return [v for v in BoundingBox if ParallelepipedContains(v)]



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
    def __init__(s, M, v0=None, use_coxiter=False):
        s.use_coxiter = use_coxiter
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
        #s.mu = (sorted(s.M1.eigenvalues()))[0]
        
        s.W=[s.V(w) for w in GetIntegerPoints(s.V1.matrix().insert_row(s.n-1, s.v0))] # w_1..w_m
        s.W.sort(key = lambda x: -s.v0.inner_product(x))

        s.En=abs(s.M.det()/GCD_list(s.M.adjoint().list()))
        s.root_lengths = [k for k in range(1,2*s.En+1) if ((2*s.En)%k == 0)]

        #print([v0.inner_product(w) for w in W])
        s.check_validity()
        s.roots = []
        print("Vinberg algorithm initialized for matrix \n{}\n".format(M))
        print(latex(M))


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
            if s.use_coxiter == True:
                M = [[ t.inner_product(r) for t in s.roots] for r in s.roots]
                print('checking polyhedron with Gram matrix')
                print(matrix(M))
                return coxiter.run(M, s.n)
            polycone=[matrix(s.M.inverse()*a) for a in Cone(s.roots).dual().rays()]
            return all(q*s.M*q.transpose()<=0 for q in polycone)
    
    #@timeit
    def FindRoots(s):
        s.roots = s.FundCone()
        for root in s.NextRoot():
            s.roots.append(root)
            #print('roots found: {0}, they are:\n{1}'.format(len(s.roots),s.roots))
            if s.is_FundPoly():
                s.print_roots()
                return

    def print_roots(s):
        polycone=[matrix(s.M.inverse()*a) for a in Cone(s.roots).dual().rays()]
        vid = len([q for q in polycone if q*s.M*q.transpose()==0])
        #print('vid =', vid)
        v_fin = len([q for q in polycone if q*s.M*q.transpose()<0])
        o=len(s.roots)
        print('V_fin =', v_fin)
        print('       ')
        print('Fundamental Coxeter polyhedron with', o, 'facets,', v_fin, 'finite and', vid, 'ideal vertices is constructed; and its roots are:')
        print('       ')
        
        for z in range(o):
            print("e",end='')
            print(z+1,'=',s.roots[z])
        #print(s.roots)
        print('       ')
        print('Ideal vertices are:')
        print('       ')
        for p in range(vid):
            print([q for q in polycone if q*s.M*q.transpose()==0][p])
        print('       ')
        print(latex(s.roots))


    def FundCone(s):
        V1_roots = [v for k in s.root_lengths for v in s.Roots_decomposed_into(s.V([0]*s.n), k) if s.IsRoot(v)]
        print('roots in V1: {}'.format(V1_roots))
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
        print('FundCone returned',[s.V(r) for r in cone.dual().rays()])
        print(latex([s.V(r) for r in cone.dual().rays()]))
        return [s.V(r) for r in cone.dual().rays()]
    
        
        
    def Roots_decomposed_into(s, a, k): #k is desired inner square, a is a non-V1 component
        '''
        Here we solve the equation (a+v1, a+v1) == k for a vector v1 in V1.
        We take a vector x in the basis g of V1, so v1 = g.x, and expand the equation to the form
        ( 2 a M g  + x^t g^t M g ) x == k - a M a^t,
        or, introducing m1, m2, c:
        ( 2 m1 + x^t m2 ) x == c.
        '''
        g = Matrix(s.V1.gens()) 
        m2 = np.dot( np.dot(g, M), g.transpose())
        m1 = np.dot( np.dot( Matrix(a), M), g.transpose() )
        c = int(k - a.inner_product(a))
        solutions = qsolve.qsolve(m2, int(2)*m1, -c)
        def V1_vector(x):
            return sum(x[i]*s.V1.gens()[i] for i in range(s.n-1))
        return [V1_vector(x)+a for x in solutions]

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
          #print(a, k, -a.inner_product(s.v0)/math.sqrt(k))
        new_roots = [v for v in s.Roots_decomposed_into(a, k) if s.IsNewRoot(v)]
        if len(new_roots)>0:
            print('new root candidates', new_roots)
        for root in new_roots:
            yield root



if __name__ == "__main__":
    #t0 = time.time()            

    # M is an inner product (quadratic form), v0 is a chosen vector

    U = matrix([[0,1],[1,0]]) # a standard hyperbolic 2-dim lattice

    D4 = matrix([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]]) # D4 lattice


    # M = diagonal_matrix(ZZ,[-15,1,1,1])
    # M = matrix([[-1,0,0,0],[0,2,-1,0],[0,-1,2,-1],[0,0,-1,2]])
    # M = matrix([[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]])
    # M = diagonal_matrix(ZZ,[-3,5,1,1])
    # M = diagonal_matrix(ZZ,[-1,3,3,2])
    # M = block_diagonal_matrix(matrix([-1]),D4,D4,D4)
    # M = matrix([[-30,0,0,0],[0,1,0,0],[0,0,2,1],[0,0,1,2]])

    # M = block_diagonal_matrix(matrix([-1]),D4)

    M = block_diagonal_matrix(U,D4)

    print('initializing a VinAl instance at a variable "A"\n')
    A = VinAl(M)
    A.FindRoots()
    #print('time_final =', time.time() - t0)

