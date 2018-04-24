import math
import _equations
import numpy as np
from timeit import timeit

#print(_equations.dE(1,1,1,1,[1,2],1,1,[1],[0.25,0.25],[0.5],2,1))
_equations.initialize()
a=np.array(np.array([1,2.,3]))
print(_equations.dTest(a))
#print(_equations.np_exp(a))
print(_equations.np_dot(a,a))
#print(_equations.np_sum(a))
quit()
