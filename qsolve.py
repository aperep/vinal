#import numba
import numpy as np
import itertools
import time
import math
import scipy.sparse.linalg as linalg

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
'''

def qform_minimum(A, m1):
  b=-.5*m1
  x = linalg.cg(A,b)[0]
  return x, np.dot(x,np.dot(A,x))+np.dot(m1,x)


def qsolve_iterative(m2, m1, c):
    n = len(m2)
    if n==1:
        a=m2[0,0]
        b=m1[0]
        disc = b*b-4*a*c
        if disc<0:
            return None
        if disc==0:
            x = int((-b)//(2*a))
            return [[x,],] if a*x*x+b*x+c==0 else []
        disc = round(math.sqrt(disc))
        x1 = int((-b+disc)//(2*a))
        x2 = int((-b-disc)//(2*a))
        sols = [[x,] for x in (x1,x2) if a*x*x+b*x+c==0]
        return sols

    x, val = qform_minimum(m2,m1)
    #print(m2, m1, c, x, val)
    if val>-c+.1:
      return None
    
    def M2(a):
        return m2[:n-1,:n-1]
    def M1(a):
        return m1[:n-1]+a*m2[n-1,:n-1]+a*m2[:n-1,n-1]
    def C(a):
        return m2[n-1,n-1]*a*a+m1[n-1]*a + c
    def oneside(a, direction):
        a=int(a)
        sols=[]
        sols_a = []
        while sols_a!=None:
          #print 'trying a=',a
          sols_a = qsolve_iterative(M2(a), M1(a), C(a))
          if sols_a==None:
            break
          sols+=[s + [a,] for s in sols_a]
          #print 'solutions', sols_a
          a+=direction
        return sols
    
    sols=oneside(math.floor(x[n-1]),-1)+oneside(math.floor(x[n-1])+1,1)
    return sols

#@timeit
def qsolve(m2, m1, c):
    s = qsolve_iterative(np.array(m2), np.array(m1).reshape(-1), c)
    return s



if __name__ == "__main__":
    tests = [
             ([[1,0,0],[0,1,0],[0,0,1]], [0,0,0], 9, 10),
             ([[3,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], [0,0,0,0], -9, 10),
    ]
    for m2, m1, c, b in tests:
        print(timeit(qsolve)(m2, m1, c))
    qform_minimum(np.array([[ 2, -1],[-1, 2]]), np.array([ 0, 128]))
