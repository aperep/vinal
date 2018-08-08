from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("qsolve.pyx", annotate=True, compiler_directives={'linetrace': True}),
)