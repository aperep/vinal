import sympy
import cdd # there is a problem with installing 'pip3 install pycddlib', since cddlib uses GMP in a broken way, see issue https://github.com/mcmtroffaes/pycddlib/issues/2 ; may try to install GMP via https://www.mersenneforum.org/showthread.php?t=23079

'''
using pycddlib for double description of polyhedra;
for self-contained implementation, look at page 5-8 of http://www.lab2.kuis.kyoto-u.ac.jp/~avis/courses/pc/2010/notes/lec5.pdf or function `initial_pair` at https://github.com/sagemath/sage/blob/e8633b09919542a65e7e990c8369fee30c7edefd/src/sage/geometry/polyhedron/double_description.py
for algorithm, see ftp://ftp.math.tu-berlin.de/pub/Lehre/LinOpt/WS06/Uebung15.11.06/fukuda96double.pdf
'''

def dot(a,b):
    return sum(a[i]*b[i] for i in range(len(a)))

def dual(rays):
    A = [[0]+list(r) for r in rays]
    A = cdd.Matrix(A)
    A.rep_type = cdd.RepType.INEQUALITY
    P = cdd.Polyhedron( cdd.Matrix(A) )
    generators = P.get_generators()
    #print(generators, '\n ___ \n')    
    dual_rays = [r[1:] for r in generators if r[0]==0]
    for i in generators.lin_set:
        dual_rays.append(tuple(-c for c in generators[i][1:]))
    return dual_rays

class Cone:
    def __init__(self, rays):
        if len(rays)==0: # for an empty cone, use the zero ray input: [[0]*n]
              print("Cone ERROR: empty input\n")
              raise ValueError
        self.rays = {tuple(r) for r in rays}
        # don't self.reduce() in init, infinite loop occurs

    def append(self, ray):
        self.rays.add(tuple(ray))
        self.reduce()
        
    def reduce(self):
        while self.reduce_once() == True:
            pass

    def reduce_once(self):
        if len(self.rays)<2: # TODO: check for a one-element list with a zero ray
          return False
        for r in self.rays:
            rays_except_r = self.rays.difference({r})
            if Cone(rays_except_r).is_nonnegative(r):
                self.rays.remove(r)
                return True
        else:
          return False

    def is_nonnegative(self, ray):
        dual_rays = dual(self.rays)
        #print('negativity test, rays:\n',self.rays, "\ndual:\n", dual_rays, '\ntest ray\n', ray)
        #print('result:', all(dot(p_ray,ray)>=0 for p_ray in dual_rays))
        return all(dot(p_ray,ray)>=0 for p_ray in dual_rays)
            
    def intersects(self, ray):
        dual_rays = dual(self.rays)
        return not( all(dot(p_ray,ray)>=0 for p_ray in dual_rays) or all(dot(p_ray,ray)<=0 for p_ray in dual_rays) )

    def __repr__(self):
        return self.rays.__repr__()


if __name__ == '__main__':
    A = [
        [1,1,1],
        [1,1,-1],
        [1,-1,1],
        [1,-1,-1],
        [1,0,0],
    ]
    B = [
        [0,1,0],
    ]
    print('Dual to [0,1,0]:',dual(B))
    C = Cone(A)
    C.reduce()
    print(C)
    print(dual([[0,0,0]]))
    print('End')