from setuptools import setup
import numpy
from Cython.Build import cythonize

setup(
    name='Integrate in 3d',
    ext_modules = cythonize('integrate_3d.pyx',annotate=True,
                            compiler_directives={'boundscheck' : False,
                                                 'wraparound' : False,
                                                 'nonecheck' : False,
                                                 'cdivision' : True}),
    include_dirs=[numpy.get_include()],
    zip_safe=False,
)
