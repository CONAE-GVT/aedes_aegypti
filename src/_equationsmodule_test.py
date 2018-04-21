import math
import _equations
from timeit import timeit

print(_equations.dE(1,1,1,1,[1,2],1,1,[1],[0.25,0.25],[0.5],2,1))
quit()

print('performance')
N=2279280
print('c++ object, R_D:  %s'%timeit('_equations.getInstance().R_D(0,299.15)', number=N, setup='import _equations'))
print('c++, R_D: %s'%timeit('R_D(0,299.15)', number=N, setup='from _equations import R_D'))
print('python, R_D: %s'%timeit('R_D(0,299.15)', number=N, setup='from __main__ import R_D'))
