from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy as np
import sys

sys.path.append('../')

setup(
    name = 'nu_funcs',
    ext_modules = cythonize([Extension('s_and_q_funcs',
                                       ['s_and_q_funcs.pyx',
                                        '../s1.c','../s2.c',
                                        '../q1.c','../q2.c',
                                        '../chbevl.c'])],
                            compiler_directives =
                            {'boundscheck' : False}),
    include_dirs = [np.get_include(),'../'],
)
