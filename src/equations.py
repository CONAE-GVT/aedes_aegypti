#coding: utf-8
import numpy as np
import math

vR_D_298K=np.array([0.24,0.2088,0.384,0.216,0.372])
#ro_25_=[0.01066,0.00873,0.01610,0.00898] #replaced by R_D_298K. which:  R_D_298K ~ ro_25*24  #(24 hours)
vDeltaH_A=np.array([10798.0,26018.0,14931.0,15725.0,15725.0])
vDeltaH_H=np.array([100000.0,55990.0,-472379.00,1756481.0,1756481.0]) #-472379 vs. -473379
vT_1_2H=np.array([14184.0,304.6,148.0,447.2,447.2])


#<precipitation related functionality v>

#Ivanov
def QR(RH_t,T_t):#in cm/day
    return 6e-5*(25 + T_t-273.15)**2 * (100.-RH_t) * 0.1#cm

def QG(p_t):#Quantity gathered#in cm
    return p_t*0.1#cm

def dW(vBS_h,vW,p_t,RH_t,T_t):#in cm/day
    netQ=QG(p_t)-QR(RH_t,T_t)
    return np.minimum(vBS_h-vW,np.maximum(netQ,-vW))

def a0(W):
    return 70.0* W

def vGamma(vL,vBS_a,vW):
    return np.array([gamma(vL[i],vBS_a[i],vW[i]) for i in range(0,len(vL))])

#</precipitation related functionality v>


def vR_D(T_t):#day^-1
    R=1.987 # universal gas constant
    return vR_D_298K * (T_t/298.0) * np.exp( (vDeltaH_A/R)* ((1.0/298.0)- (1.0/T_t)) ) / ( 1.0+ np.exp( (vDeltaH_H/R)* ((1.0/vT_1_2H)-(1.0/T_t)) ) )


def gamma(L,BS,W):
    epsilon=1e-4
    if(BS==0 or W <epsilon):#W *1000./BS_s <0.1
        return 1.0#no water total inhibition#1960 Aedes aegypti (L.), The Yellow Fever Mosquito(Page 165)
    if(L/BS<=a0(W)-epsilon):
        return 0
    elif(a0(W)-epsilon < L/BS <a0(W)+epsilon):
        #a (a0-e) + b=0 => b=-a (a0 -e)
        #a (a0 + e) + b=0.63 => a(a0+e) - a(a0-e) = 2 a e = 0.63 =>a=0.63/(2 e)
        a=0.63/(2.0*epsilon)
        b=-a*(a0(W)-epsilon)
        return a * (L/BS) + b
    elif(L/BS>=a0(W)+epsilon):
        return 0.63


def f(vW,vBS_d):#TODO:change name to something with meaning
    epsilon=1e-4
    vf=vW/(vW+epsilon) * vBS_d
    if(vf.max()<1e-20):
        return vf
    else:
        return vf/np.sum(vf)#TODO: check this

def dvE(vE,vL,A1,A2,vW,T_t,BS_a,vBS_d,elr,ovr1,ovr2):
    egn=63.0
    me=0.01#mortality of the egg, for T in [278,303]
    return egn*( ovr1 *A1  + ovr2* A2)*f(vW,vBS_d) - me * vE - elr* (1-vGamma(vL,BS_a*vBS_d,vW)) * vE

def dvL(vE,vL,vW,T_t,BS_a,vBS_d,elr,lpr,vAlpha0):
    ml=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#mortality of the larvae, for T in [278,303]
    vAlpha=vAlpha0/(BS_a*vBS_d)
    return elr* (1-vGamma(vL,BS_a*vBS_d,vW)) * vE - ml*vL - vAlpha* vL*vL - lpr *vL #-35.6464*(1.-beta(vW,vBS_od,vBS_id))*L#-24.*(1.-beta(vW))*L# -log(1e-4/5502.)/(1.)=17.823207313460703

def dvP(vL,vP,T_t,lpr,par):
    mp=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#death of pupae
    return lpr*vL - mp*vP  - par*vP

def dA1(vP,A1,T_t,par,ovr1):
    ef=0.83#emergence factor
    ma=0.09#for T in [278,303]
    return np.sum(par*ef*vP/2.0) - ma*A1 - ovr1*A1

def dA2(A1,A2,T_t,ovr1):
    ma=0.09#for T in [278,303]
    return ovr1*A1 - ma*A2

def diff_eqs(Y,t,parameters):
    '''The main set of equations'''
    T_t=parameters.weather.T(t)
    p_t=parameters.weather.p(t)
    RH_t=parameters.weather.RH(t)
    elr,lpr,par,ovr1,ovr2=vR_D(T_t)
    BS_a,vBS_h,vBS_s,vBS_d,vAlpha0,n=parameters.BS_a,parameters.vBS_h,parameters.vBS_s,parameters.vBS_d,parameters.vAlpha0,parameters.n
    EGG,LARVAE,PUPAE,ADULT1,ADULT2,WATER=parameters.EGG,parameters.LARVAE,parameters.PUPAE,parameters.ADULT1,parameters.ADULT2,parameters.WATER

    vE,vL,vP,A1,A2,vW=Y[EGG],Y[LARVAE],Y[PUPAE],Y[ADULT1],Y[ADULT2],Y[WATER]

    dY=np.zeros((3*n + 2 + n ))
    dY[EGG]    = dvE(vE,vL,A1,A2,vW,T_t,BS_a,vBS_d,elr,ovr1,ovr2)
    dY[LARVAE] = dvL(vE,vL,vW,T_t,      BS_a,vBS_d,elr,lpr,vAlpha0)
    dY[PUPAE]  = dvP(vL,vP,T_t,lpr,par)
    dY[ADULT1] = dA1(vP,A1,T_t,par,ovr1)
    dY[ADULT2] = dA2(A1,A2,T_t,ovr1)
    dY[WATER]  = dW(vBS_h,vW,p_t,RH_t,T_t)

    return dY   # For odeint
