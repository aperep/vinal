import sys
sys.path.insert(0, '../src/sage')
from vinal import *

def test_U_D4():
    U = matrix([[0,1],[1,0]]) # a standard hyperbolic 2-dim lattice
    D4 = matrix([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]]) # D4 lattice
    M = block_diagonal_matrix(U,D4)
    A = VinAl(M,output=None)
    A.FindRoots()
    assert len(A.roots) == 6
    return A.roots