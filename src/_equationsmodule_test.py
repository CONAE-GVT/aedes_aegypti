import math
import _equations
from timeit import timeit

vR_D_298K=[0.24,0.2088,0.384,0.216,0.372]
#ro_25_=[0.01066,0.00873,0.01610,0.00898] #replaced by R_D_298K. which:  R_D_298K ~ ro_25*24  #(24 hours)
vDeltaH_A=[10798.0,26018.0,14931.0,15725.0,15725.0]
vDeltaH_H=[100000.0,55990.0,-472379.00,1756481.0,1756481.0] #-472379 vs. -473379
vT_1_2H=[14184.0,304.6,148.0,447.2,447.2]

def sum(a,b):
    return a+b

def mul(a,b):
    return a*b

def R_D(stage,T_t):#day^-1
    R=1.987 # universal gas constant
    R_D_298K=vR_D_298K[stage]
    #ro_25=ro_25_[stage]
    deltaH_A=vDeltaH_A[stage]
    deltaH_H=vDeltaH_H[stage]
    T_1_2H=vT_1_2H[stage] # K
    return R_D_298K * (T_t/298.0) * math.exp( (deltaH_A/R)* ((1.0/298.0)- (1.0/T_t)) ) / ( 1.0+ math.exp( (deltaH_H/R)* ((1.0/T_1_2H)-(1.0/T_t)) ) )



equations = _equations.getInstance()
#print(equations.sum(3,5))
print(equations.mul(3,5))#from an object
print(_equations.mul(3,5))#static
print(mul(3,5))

print('performance')
N=2279280
print('c++ object, mul: %s'%timeit('_equations.getInstance().mul(3,5)', number=N, setup='import _equations'))
print('c++, mul: %s'%timeit('mul(3,5)', number=N, setup='from _equations import mul'))#it seems to be faster than "import _equations; _equations.mul"
print('python, mul: %s'%timeit('mul(3,5)', number=N, setup='from __main__ import mul'))

print('c++ object, R_D:  %s'%timeit('_equations.getInstance().R_D(0,299.15)', number=N, setup='import _equations'))
print('c++, R_D: %s'%timeit('R_D(0,299.15)', number=N, setup='from _equations import R_D'))
print('python, R_D: %s'%timeit('R_D(0,299.15)', number=N, setup='from __main__ import R_D'))
