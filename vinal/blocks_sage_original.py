from sage.all import *
import qsolve_sage as qsolve
import sympy as sp

def fundamental_cone(s, diag=False):
        V1_roots = [vector([int(i) for i in list(v.T)]) for k in s.root_lengths for v in s.roots_of_type(zero_vector(s.n), k) if s.is_root(v)]
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
        return sp.Matrix([r for r in cone.dual().rays()]).T

                
def roots_of_type(s, a, k, diag=False): #k is desired inner square, a is a non-V1 component
        '''
        Here we solve the equation (a+v1, a+v1) == k for a vector v1 in V1.
        We take a vector x in the basis g of V1, so v1 = g.x, and expand the equation to the form
        ( 2 a M g  + x^t g^t M g ) x == k - a M a^t,
        or, introducing m1, m2, c:
        ( 2 m1 + x^t m2 ) x == c.
        '''
        
        if diag == True:
          raise NotImplementedError
        g = s.V1_basis_diag.T
        print(f'shapes g {g.shape} s.Q {s.Q.shape}')
        m2 = g * s.Q * g.transpose()
        m1 = sp.Matrix(a).T * s.Q * g.transpose()
        c = int(k - a.inner_product(a))
        solutions = qsolve.qsolve(m2, int(2)*m1, -c)
        print('solutions',solutions)
        def V1_vector(x):
            return sp.Matrix(x).T * g
        return [V1_vector(x).T+sp.Matrix(a) for x in solutions]

        
def root_types(s, stop=-1): # iterates pairs (w_i + c v_0, ||a||) from minimum, infinity or `stop` times
        candidates = {k:1 for k in s.root_lengths} # dictionary of vector numbers for each k: we have a series of vectors {W}\0, {v0 + W}, etc. The number is the position in this series.
        def cand_a(n): # the non-V1 component of vector under number n
            return s.W[n%len(s.W)]+ (n//len(s.W))*s.v0
        while True:
            k = min(s.root_lengths, key=lambda k: -s.v0.inner_product(cand_a(candidates[k]))/math.sqrt(k)) # can be optimized with heapq
            yield cand_a(candidates[k]), k
            candidates[k]+=1
            stop-=1
            if stop == 0:
                return
