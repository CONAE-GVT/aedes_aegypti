#>python setup.py build_ext --inplace
#python -m timeit -s 'from equations import R_D' 'R_D(1,315.15)'
#python -m timeit -s 'from _equations import R_D' 'R_D(1,315.15)'
from distutils.core import setup
from Cython.Build import cythonize
import os
import shutil
os.environ['CFLAGS'] = '-O3'
shutil.copyfile('equations.py', '_equations.pyx')
setup(
    ext_modules=cythonize('_equations.pyx'),
)
