import numpy as np
import _equations
_equations.initialize()
from _equations import vR_D,dvE,dvL,dvP,dA1,dA2,dW

EGG=0
LARVAE=1
PUPAE=2
ADULT1=3
ADULT2=4
WATER=5

def diff_eqs(Y,t,parameters):
    '''The main set of equations'''
    T_t=float(parameters.weather.T(t))
    p_t=parameters.weather.p(t)
    ws_t=parameters.weather.ws(t)
    wss_t=parameters.ws_s*ws_t
    elr,lpr,par,ovr1,ovr2=[float(r) for r in vR_D(T_t)]
    BS_a,vBS_oc,vBS_ic,vBS_d,vBS_os,n,m=parameters.BS_a,parameters.vBS_oc,parameters.vBS_ic,parameters.vBS_d,parameters.vBS_os,parameters.n,parameters.m
    EGG,LARVAE,PUPAE,ADULT1,ADULT2,WATER=range(0,n+m),range(n+m,2*(n+m)),range(2*(n+m),3*(n+m)),3*(n+m),3*(n+m)+1,3*(n+m)+2

    vE,vL,vP,A1,A2,vW=Y[EGG],Y[LARVAE],Y[PUPAE],Y[ADULT1],Y[ADULT2],Y[WATER:]
    vW=np.concatenate((vW,vBS_ic))#constant water for inside BS

    dY=np.zeros((3*(n+m)+2+n))
    dY[EGG]    = dvE(vE,vL,A1,A2,vW,T_t,BS_a,vBS_d,elr,ovr1,ovr2)
    dY[LARVAE] = dvL(vE,vL,vW,T_t,      BS_a,vBS_d,elr,lpr)
    dY[PUPAE]  = dvP(vL,vP,T_t,lpr,par)
    dY[ADULT1] = dA1(vP,A1,T_t,par,ovr1)
    dY[ADULT2] = dA2(A1,A2,T_t,ovr1)
    dY[WATER:] = [dW(vW[i],vBS_oc[i],vBS_os[i],T_t,p_t,wss_t,t) for i in range(0,n)]#note that this goes till n.

    return dY   # For odeint
