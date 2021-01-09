from vinal import *

D4 = Matrix([[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]])
U = Matrix([[0,1],[1,0]])
V = VinAl(diag(U,D4))
roots = [tuple(r) for r in V.roots_in_v0_perp()]
v = (0, 0, -1, -1, 0, 0)
print('roots in V1 are', roots)
print(f'roots in V1 should include e.g. v = {v}')
v = Matrix(v)
print(f'v is root of length 2 orthogonal to v0: V.is_root(v)={V.is_root(v)}, length={V.dot(v,v)}, (v,v0)={V.dot(v,V.v0)}')

print(f'diagonal basis vectors should be primitive: {V.basis_diag.T}')