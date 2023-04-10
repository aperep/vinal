import time
from line_profiler import LineProfiler
import cProfile

def timeit(f): # @timeit 
    
    def timed(*args, **kw):
        
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        
        print 'func:%r args:[%r, %r] took: %2.4f sec' % \
            (f.__name__, args, kw, te-ts)
        return result    
    return timed


# taken from https://zapier.com/engineering/profiling-python-boss/

def timefunc(f): # @timefunc 
    def f_timer(*args, **kwargs):
        start = time.time()
        result = f(*args, **kwargs)
        end = time.time()
        print f.__name__, 'took', end - start, 'time'
        return result
    return f_timer

class timewith(): # with timewith('%chunk_name%') as timer: %chunk% timer.checkpoint('done with %chunk%')
    def __init__(self, name=''):
        self.name = name
        self.start = time.time()

    @property
    def elapsed(self):
        return time.time() - self.start

    def checkpoint(self, name=''):
        print '{timer} {checkpoint} took {elapsed} seconds'.format(
            timer=self.name,
            checkpoint=name,
            elapsed=self.elapsed,
        ).strip()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.checkpoint('finished')
        pass


def do_cprofile(func): # @do_cprofile
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func


def do_profile(follow=[]): # @do_profile(follow=[%inner_func%])
        def inner(func):
            def profiled_func(*args, **kwargs):
                try:
                    profiler = LineProfiler()
                    profiler.add_function(func)
                    for f in follow:
                        profiler.add_function(f)
                    profiler.enable_by_count()
                    return func(*args, **kwargs)
                finally:
                    profiler.print_stats()
            return profiled_func
        return inner