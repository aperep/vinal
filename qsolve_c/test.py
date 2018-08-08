from subprocess import call, check_output
import time
import numpy as np
import warnings
import cProfile

import line_profiler

def timeit(f): # @timeit
    
    def timed(*args, **kw):
        
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        
        print 'func:%r args:[%r, %r] took: %2.4f sec' % \
            (f.__name__, args, kw, te-ts)
        return result
    return timed

    
if __name__ == "__main__":
    warnings.filterwarnings("ignore", message="numpy.dtype size changed")
    call('bash -c "python setup.py build_ext --inplace"', shell=True)
    import qsolve
    tests = [
             ([[1,0,0],[0,1,0],[0,0,1]], [0,0,0], 9, 10),
             ([[3,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], [0,0,0,0], -9, 10),
             ([[3,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]], [0,0,0,0,0], -9, 10),
             ([[3,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,1,1]], [0,0,0,0,0,1], -9, 15),
    ]
    l = line_profiler.LineProfiler()
    l.add_function(qsolve.qsolve)
    for m2, m1, c, b in tests:
        timeit(qsolve.qsolve)(m2, m1, c)
        cProfile.run('qsolve.qsolve(m2, m1, c)', sort='time')
        l.run('qsolve.qsolve(m2, m1, c)')
        l.print_stats()

        #qsolve.qform_minimum(np.array([[ 2, -1],[-1, 2]]), np.array([ 0, 128]))
