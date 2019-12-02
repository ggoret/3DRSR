
from distutils.core import setup
from Cython.Build import cythonize

openmp_option = ["-fopenmp"]

setup(
    name = 'MDSP Cython module',
    ext_modules = cythonize("rotsum.pyx", extra_compile_args = openmp_option, extra_links_args = openmp_option, annotate=True),
)
