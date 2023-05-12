import sys
sys.path.insert(0, '../src/sage')
from vinal import *

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
