import sys
import datetime
import numpy as np
from scipy import interpolate
import otero_precipitation as op

current_module = sys.modules[__name__]
Wip=None
def waterEquations(vW,t):
    T_t=op.T(t)
    vW[vW<0]=0#this is to make rk work
    dW = [op.dW(vW[i],op.vBS_oc[i],op.vBS_os[i],T_t,t) for i in range(0,op.n)]
    return dW

def populationEquations(Y,t):
    '''The main set of equations'''
    dY=np.zeros((5))
    Y[Y<0]=0#this is to make rk work
    T_t=op.T(t)
    E,L,P,A1,A2=Y
    vW=np.array(Wip(t))
    dY[op.EGG]    = op.dE(E,L,A1,A2,vW,T_t)
    dY[op.LARVAE] = op.dL(E,L,vW,T_t)
    dY[op.PUPAE]  = op.dP(L,P,T_t)
    dY[op.ADULT1] = op.dA1(P,A1,T_t)
    dY[op.ADULT2] = op.dA2(A1,A2,T_t)
    return dY   # For odeint


def solveWaterEquations(initial_condition =[0. for i in range(0,op.n)]):
    #print('Water Equations start:%s')%datetime.datetime.now().time()
    time_range,W0,W=op.solveEquations(initial_condition=initial_condition,equations=waterEquations,method='odeint')
    #print('Water Equations end:  %s')%datetime.datetime.now().time()
    vip=[interpolate.InterpolatedUnivariateSpline(time_range,W[:,i]) for i in range(0,op.n)]
    current_module.Wip=lambda t:[vip[i](t) for i in range(0,op.n)]
    return time_range,W0,W

def solvePopulationEquations(initial_condition = [100.0, 0.0,0.0,0.0,0.0]):
    #print('Population Equations start:%s')%datetime.datetime.now().time()
    time_range,Y0,Y=op.solveEquations(initial_condition=initial_condition,equations=populationEquations,method='rk')
    #print('Population Equations end:  %s')%datetime.datetime.now().time()
    return time_range,Y0,Y

def solveEquations(initial_condition = [100.0, 0.0,0.0,0.0,0.0]+ [0. for i in range(0,op.n)]):
    time_range,W0,W=solveWaterEquations(initial_condition = initial_condition[op.WATER:op.WATER+op.n])
    time_range,Y0,Y=solvePopulationEquations(initial_condition=initial_condition[0:op.WATER])

    #backward compatibility
    RES=np.zeros((len(time_range),5+op.n))
    RES[:,:op.WATER]=Y
    RES[:,op.WATER:]=W

    return time_range,initial_condition,RES
