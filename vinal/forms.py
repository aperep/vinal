from sympy import *
import itertools
import math
import functools
from fractions import Fraction

def negative_vector(Q):
    D, T = rational_diagonal_form(Q)
    #print('\n'.join(["M=",repr(M),"D=",repr(D),"T=",repr(T)]))
    #print [(D[i][i], T.column(i)) for i in range(D.nrows())]
    vectors = [T.col(i) for i in range(D.shape[0]) if D[i,i]<0]
    assert len(vectors) == 1
    v0 = vectors[0]
    v0 = v0 / gcd(*v0)
    return Matrix(v0)


def rational_diagonal_form(Q, d=None):
    """
    takes a quadratic form Q (of type sympy.Matrix) and returns its diagonalization D in coordinates T (both in integers). If Q is hyperbolic, it puts negative vector on the first place.
    This is an adapted version of the eponymous function from the SageMath project https://github.com/sagemath/sagelib/blob/master/sage/quadratic_forms/quadratic_form__local_field_invariants.py
    """
    n = Q.shape[0]
    D = Q[:,:]


    # make symmetric
    for i in range(n):
      for j in range(i+1,n):
        D[i,j] = D[j,i] = (D[j,i] + D[i,j])/2
    
    T = eye(n)
    init_printing()
    ## Clear the entries one row at a time.
    for i in range(n):
      

        ## Deal with rows where the diagonal entry is zero.
        if D[i,i] == 0:

            ## Look for a non-zero entry and use it to make the diagonal non-zero (if it exists)
            for j in range(i+1, n):
                if D[i,j] != 0:
                    temp = eye(n)
                    if D[i,j] + D[j,j] == 0:
                        temp[j, i] = -1
                    else:
                        temp[j, i] = 1

                    ## Apply the transformation
                    
                    D = temp.T*D*temp
                    T = T * temp
                    break

        ## Create a matrix which deals with off-diagonal entries (all at once for each row)
        temp = eye(n)
        for j in range(i+1, n):
            if D[i,j] != 0: ## This should only occur when Q[i,i] != 0, which the above step guarantees.
                temp[i,j] = -D[i,j]    
                temp[j,j] =  D[i,i] 
        D = temp.T*D*temp
        T = T * temp

    # put negative vector first
    if any(D[i,i]<0 for i in range(n)):    
      negative_index = min(i for i in range(n) if D[i,i]<0)
      if negative_index != 0:
        temp = eye(n)
        temp[0, 0], temp[negative_index, negative_index] = 0, 0
        temp[0, negative_index], temp[negative_index, 0] = 1, 1
        D = temp.T*D*temp
        T = T * temp

    # make basis vectors primitive:
    for i in range(n):        
        GCD = functools.reduce(lambda x,y:gcd(x,y),tuple(T[:,i]))
        if abs(GCD) != 1:
            temp = eye(n)
            temp[i,i]=Fraction(1,GCD)
            D = temp.T*D*temp
            T = T * temp

    return simplify(D), simplify(T)



if __name__ == '__main__':
  D = diag(-3,5,1,1)
  Q = rational_diagonal_form(D)
  print(Q)