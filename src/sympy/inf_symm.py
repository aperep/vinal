import itertools
from itertools import product
import numpy as np

M = np.diagonal_matrix(ZZ,[-23,1,1,1])

roots = [[0, -1, -1, 0], (0, 0, 1, -1), (0, 1, 0, 0), (1, 0, 3, 4), (1, 0, 0, 5), (1, -2, 2, 4), (2, -3, 6, 7), (2, -2, 3, 9), (3, -8, 8, 9), (4, -9, 12, 12), (6, -1, 10, 27), (7, -2, 10, 32), (10, -22, 27, 33)]

print (np.matrix(roots[0]))

#B = product(*[g for i in range(4)])

#D = product(*[g for i in range(4)])

print('Equiv Matrices:')
'''
for b1 in B:
    X1 = [[inner(u,v).real() for u in b1] for v in b1]
    if (np.linalg.det(np.matrix(X1)) != 0):
        print(X1)
'''

def inner(u,v):
  return (np.matrix(u)*np.matrix(M)*np.matrix(v).transpose())

print (inner(roots[0],roots[1]))


'''
for b1 in B:
    X1 = [[inner(u,v).real() for u in b1] for v in b1]
    if (np.linalg.det(np.matrix(X1)) != 0):
        #print(X1)
        for b2 in D:
            #X1 = [[inner(u,v).real() for u in b1] for v in b1]
            X2 = [[inner(u,v).real() for u in b2] for v in b2]
            #print(X1)
            #print(np.linalg.det(np.matrix(X1)))
            if(X1 == X2):
                print('000000000000000000')
                print('b1 =', b1)
                print('b2 =', b2)
                Z = np.matrix([[Mult(b2,Inv(b1))[i][j].real() for i in range(4)] for j in range(4)])
                print(np.matrix(Mult(b2,Inv(b1))))
                print(Det(Mult(b2,Inv(b1))))
                w,v = np.linalg.eig(Z)
                print('w =', w)
                print('v =', v)
                print('Mult2 =', Mult(Mult(b2,Inv(b1)),Mult(b2,Inv(b1))))

'''

