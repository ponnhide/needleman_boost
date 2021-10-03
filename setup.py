from distutils.core import setup, Extension
from Cython.Build import cythonize
from numpy import get_include 

ext = Extension("_global_alignment", sources=["_global_alignment.pyx"], include_dirs=['.', get_include()], extra_compile_args=["-O3"], language="c++")
setup(name="_global_alignment", ext_modules=cythonize([ext]))
