import numpy as np
from _equations import dE,dL,dP,dA1,dA2,dW

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
    BS_a,vBS_oc,vBS_ic,vBS_od,vBS_id,vBS_os,n,m=parameters.BS_a,parameters.vBS_oc.tolist(),parameters.vBS_ic.tolist(),parameters.vBS_od.tolist(),parameters.vBS_id.tolist(),parameters.vBS_os.tolist(),parameters.n,parameters.m
    Y=Y.tolist()
    E,L,P,A1,A2=Y[:WATER]
    vW=np.array(Y[WATER:]).tolist()

    dY=np.zeros((5+n))
    dY[EGG]    = dE(E,L,A1,A2,vW,T_t,BS_a,vBS_ic,vBS_od,vBS_id,n,m)
    dY[LARVAE] = dL(E,L,vW,T_t,      BS_a,vBS_ic,vBS_od,vBS_id,n,m)
    dY[PUPAE]  = dP(L,P,T_t)
    dY[ADULT1] = dA1(P,A1,T_t)
    dY[ADULT2] = dA2(A1,A2,T_t)
    dY[WATER:] = [dW(vW[i],vBS_oc[i],vBS_os[i],T_t,p_t,wss_t,t) for i in range(0,n)]

    return dY   # For odeint
