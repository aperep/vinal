import line_profiler
from vinal import *

V = VinAl(diag(-1,1,1,1,1,1,1,1,1,1)) # unimodular lattice of signature (9,1)
V3 = VinAl(diag(-3,5,1))
profile = line_profiler.LineProfiler()
V.fundamental_cone()


