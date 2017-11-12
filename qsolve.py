#import numba
import numpy as np
import itertools

'''
This module solves a general quadratic diophantine equation in n variables x=(x_1,..,x_n) in the form 
(x^t).M2.x + M1.x + c = 0,
where M2 is a positive definite matrix.
Optional parameter 'boundary' restricts solutions to the cube |x_1|, .. |x_n| < boundary.
If 'boundary' is present, then M2 can be arbitrary.
'''

# egcd() and crt() are taken from https://gist.github.com/sirodoht/ee2abe82eca70f5b1869
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


def vectors_multiply(vectors, primes):
        n = len(vectors[0])
        return tuple(crt(primes,[v[i] for v in vectors]) for i in range(n))



def qsolve(m2, m1, c, boundary=None):
    #print("qsolve with arguments", m2, m1, c)
    if boundary == None:
        raise NotImplementedError
    n = len(m2)
    m2 = np.array(m2)
    m1 = np.array(m1)
    primes = [2,3,5] # TODO: choose dynamically w.r.t. boundary

    def diff(x):
        x1 = np.dot(x, m2)
        ans2 = np.dot(x1+m1, x) + c
        #print(type(ans2))
        return int(ans2)

    def solve_mod_p(p):
        Box = itertools.product(*[range(p) for i in range(n)])
        return [x for x in Box if diff(x)%p == 0]
    
    sols_mod_primes = [solve_mod_p(p) for p in primes]
    sol_combinations = [i for i in itertools.product(*sols_mod_primes)]
    remainders = set([ vectors_multiply(sol, primes) for sol in sol_combinations]) # need to filter repeated vectors
    modulus = np.prod(primes) # m=modulus; plan: rebalance to [-m/2..m/2], iterate y over Box/m and look for Q(my+x)=k. Filter y's somehow.

    boundary = int(round(boundary/modulus) + 1)
    sols = []
    for y in itertools.product(*[range(-boundary,boundary) for i in range(n)]):
        for x in remainders:
            v = np.array(y)*modulus+np.array(x)
            if (diff(v) == 0):
                sols.append(v)
    return sols

                #@numba.jit(nopython=True)



if __name__ == "__main__":
    test = qsolve([[1,0,0],[0,1,0],[0,0,1]], [0,0,0], 9, 10)
    for t in test:
        print(t)