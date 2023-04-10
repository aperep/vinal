from vinal import *
import sympy

lattices = [
        {
            "name": '-2 + A2 + 1 + 1 + 1',
            "matrix": Matrix([[-2,0,0,0,0,0],[0,2,-1,0,0,0],[0,-1,2,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]]),
            "current time": 38.3,
            "reference time": 12.8,            
        },
        {
            "name": 'diag(-1,3,3,2)',
            "matrix": diag(-1,3,3,2),
            "current time": 0.94,
            "reference time": 1.4,            
        },
        {
            "name": 'diag(-3,5,1,1)',
            "matrix": diag(-3,5,1,1),
            "current time": 19.1,
            "reference time": 15.2,            
        },
        {
            "name": 'diag(-7,1,1,1)',
            "matrix": diag(-7,1,1,1),
            "current time": 1.0,
            "reference time": 1.27,            
        },
        {
            "name": 'the unimodular lattice of signature (7,1)',
            "matrix":  diag(-1,1,1,1,1,1,1,1),
            "current time": 7.8,
            "reference time": 10.5,            
        },
        {
            "name": 'the unimodular lattice of signature (8,1)',
            "matrix": diag(-1,1,1,1,1,1,1,1,1),
            "current time": 15.6,
            "reference time": 19.9,            
        },
        {
            "name": 'the unimodular lattice of signature (9,1)',
            "matrix": diag(-1,1,1,1,1,1,1,1,1,1),
            "current time": 34.9,
            "reference time": 30.8,            
        },
]


if __name__ == '__main__':
  print('The list "lattices" contains lattices their timings')