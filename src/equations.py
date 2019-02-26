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

'''
    { QG(BS_s,t)-QR(BS_s,T(t))    if 0 < W < BS_h
dW= { QG(BS_s,t)                  if W <= 0.0
    { -QR(BS_s,T(t))               if W >= BS_h
Note: in the implementation we needed to add functions to make function continuous, otherwise odeint breaks
'''
def dW(W,BS_h,T_t,p_t,RH_t):#in cm/day
    epsilon=1e-1#1mm
    if(0+epsilon < W < BS_h-epsilon):
        return QG(p_t)-QR(RH_t,T_t)
    elif(W <= 0.0+epsilon):
        return QG(p_t) - QR(RH_t,T_t)*(W/epsilon)
    elif( W >= BS_h-epsilon):
        return QG(p_t)*((BS_h-W)/epsilon) - QR(RH_t,T_t)

def a0(W):
    return 70.0* W

def vGamma(vL,vBS_a,vW):
    return np.array([gamma(vL[i],vBS_a[i],vW[i]) for i in range(0,len(vL))])

def waterEquations(vW,t,parameters):
    T_t=parameters.weather.T(t)
    p_t=parameters.weather.p(t)
    RH_t=parameters.weather.RH(t)
    vmf_t=parameters.mf(t)*parameters.vBS_mf*10.# cm -> mm
    vBS_h,n=parameters.vBS_h,parameters.n
    return [dW(vW[i],vBS_h[i],T_t,p_t+vmf_t[i],RH_t) for i in range(0,n)]
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


def ovsp(vW,vBS_d,vW_l,mBS_l):#OViposition Site Preference
    epsilon=1e-4
    vf=vW/(vW+epsilon) * vBS_d#check this method is not spontaneus generation of eggs.(like inventing adults.)
    if(vf.max()>epsilon): vf/=vf.sum()
    return np.where(mBS_l==np.floor(vW_l),1,0)*vf

def wetMask(vW_l,mBS_l):
    return np.where(mBS_l<=vW_l,1,0)

def dvE(mE,vL,A1,A2,vW_t,BS_a,vBS_d,elr,ovr1,ovr2,wet_mask,vW_l,mBS_l,egnCorrector,t):
    egn=63.0
    me=0.01#mortality of the egg, for T in [278,303]
    ovsp_t=ovsp(vW_t,vBS_d,vW_l,mBS_l)
    egn_c=egnCorrector(egn,ovr1 *A1  + ovr2* A2, t)
    return egn_c*( ovr1 *A1  + ovr2* A2)*ovsp_t - me * mE - elr* (1-vGamma(vL,BS_a*vBS_d,vW_t)) * mE*wet_mask

def dvL(mE,vL,vW,T_t,BS_a,vBS_d,elr,lpr,vAlpha0,wet_mask):
    ml=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#mortality of the larvae, for T in [278,303]
    mdl=2.#mortality of dry larvae.TODO:Unjustified!
    vAlpha=vAlpha0/(BS_a*vBS_d)
    epsilon=1e-4
    return elr* (1-vGamma(vL,BS_a*vBS_d,vW)) * np.sum(mE*wet_mask,axis=0) - ml*vL - vAlpha* vL*vL - lpr *vL -mdl*(1.- vW/(vW+epsilon))*vL

def dvP(vL,vP,T_t,lpr,par):
    mp=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#death of pupae
    return lpr*vL - mp*vP  - par*vP

def dA1(vP,A1,par,ovr1):
    ef=0.83#emergence factor
    ma=0.09#for T in [278,303]
    return np.sum(par*ef*vP/2.0) - ma*A1 - ovr1*A1

def dA2(A1,A2,ovr1):
    ma=0.09#for T in [278,303]
    return ovr1*A1 - ma*A2

def diff_eqs(Y,t,parameters):
    '''The main set of equations'''
    T_t=parameters.weather.T(t)
    elr,lpr,par,ovr1,ovr2=vR_D(T_t)
    BS_a,BS_lh,vBS_d,vAlpha0,m,n,mBS_l=parameters.BS_a,parameters.BS_lh,parameters.vBS_d,parameters.vAlpha0,parameters.m,parameters.n,parameters.mBS_l
    EGG,LARVAE,PUPAE,ADULT1,ADULT2=parameters.EGG,parameters.LARVAE,parameters.PUPAE,parameters.ADULT1,parameters.ADULT2

    vW_t=parameters.vW(t)
    vE,vL,vP,A1,A2=Y[EGG].reshape((n,m)).transpose(),Y[LARVAE],Y[PUPAE],Y[ADULT1],Y[ADULT2]
    vW_l=vW_t/BS_lh
    wet_mask=wetMask(vW_l,mBS_l)

    dY=np.empty(Y.shape)
    dY[EGG]    = dvE(vE,vL,A1,A2,vW_t,BS_a,vBS_d,elr,ovr1,ovr2,wet_mask,vW_l,mBS_l,parameters.egnCorrector,t).transpose().reshape((1,m*n))
    dY[LARVAE] = dvL(vE,vL,vW_t,T_t,BS_a,vBS_d,elr,lpr,vAlpha0,wet_mask)
    dY[PUPAE]  = dvP(vL,vP,T_t,lpr,par)
    dY[ADULT1] = dA1(vP,A1,par,ovr1)
    dY[ADULT2] = dA2(A1,A2,ovr1)

    return dY   # For odeint
