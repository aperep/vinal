#import numba
import numpy as np
import itertools
import time
import math

def timeit(f): # @timeit
    
    def timed(*args, **kw):
        
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        
        print 'func:%r args:[%r, %r] took: %2.4f sec' % \
            (f.__name__, args, kw, te-ts)
        return result
    return timed

'''
This module solves a general quadratic diophantine equation in n variables x=(x_1,..,x_n) in the form 
(x^t).M2.x + M1.x + c = 0,
where M2 is a positive definite matrix.
Optional parameter 'boundary' restricts solutions to the cube |x_1|, .. |x_n| < boundary.
If 'boundary' is present, then M2 can be arbitrary.
'''

# finds x that minimizes Q = x^T.m1.x + m2.x, see http://komarix.org/ac/papers/thesis/thesis_html/node11.html
def qform_minimum(m2, m1):
  n=len(m2)
  x=np.zeros(n)
  d=np.zeros(n)
  r=np.copy(m1)
  r_prev=r
  beta=0
  alpha=0
  while r.dot(r)>0:
    d     = r+beta*d
    alpha =-np.dot(d,b)/np.dot(d,(np.dot(A,d)))
    x     = x-alpha*d
    r_prev= r
    r     = b-np.dot(A,x)
    beta  = np.dot(r,r)/np.dot(r_prev,r_prev)
  return x, np.dot(x,np.dot(m2,x))+np.dot(m1,x)

def qsolve_iterative(m2, m1, c, boundary):
    n = len(m2)
    if n==1:
        a=m2[0,0]
        b=m1[0]
        disc = b*b-4*a*c
        if disc<0:
            return 0, None
        if disc==0:
            x = int((-b)//(2*a))
            return (1, [[x,],]) if a*x*x+b*x+c==0 else (0, [])
        disc = round(math.sqrt(disc))
        x1 = int((-b+disc)//(2*a))
        x2 = int((-b-disc)//(2*a))
        sols = [[x,] for x in (x1,x2) if a*x*x+b*x+c==0]
        return len(sols), sols

    x, val = qform_minimum(m2,m1)
    print(m2, m1, c, x, val)
    if val>-c+.1:
      return 0, None
    
    def M2(a):
        return m2[:n-1,:n-1]
    def M1(a):
        #print(m1[:n-1].shape)
        #print(m2[n-1,:n-1].shape)
        #print(m2[:n-1,n-1].shape)
        return m1[:n-1]+a*m2[n-1,:n-1]+a*m2[:n-1,n-1]
    def C(a):
        return m2[n-1,n-1]*a*a+m1[n-1]*a + c
    def oneside(a, direction):
        sols=[]
        sols_a = []
        while sols_a!=None:
          sols+=[s + [a,] for s in sols_a]
          solnum, sols_a = qsolve_iterative(M2(a), M1(a), C(a), boundary)
          a+=direction
        return sols
    
    sols=oneside(math.floor(x[n-1]),-1)+oneside(math.floor(x[n-1])+1,1)
    return len(sols), sols

#@timeit
def qsolve(m2, m1, c, boundary):
    s = qsolve_iterative(np.array(m2), np.array(m1).reshape(-1), c, boundary)[1]
    #print(s)
    return s



if __name__ == "__main__":
    tests = [
             ([[1,0,0],[0,1,0],[0,0,1]], [0,0,0], 9, 10),
             ([[3,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], [0,0,0,0], -9, 10),
    ]
    for m2, m1, c, b in tests:
        print(qsolve(m2, m1, c, b))
