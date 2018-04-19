#coding: utf-8
import numpy as np
import math

vR_D_298K=[0.24,0.2088,0.384,0.216,0.372]
#ro_25_=[0.01066,0.00873,0.01610,0.00898] #replaced by R_D_298K. which:  R_D_298K ~ ro_25*24  #(24 hours)
vDeltaH_A=[10798.0,26018.0,14931.0,15725.0,15725.0]
vDeltaH_H=[100000.0,55990.0,-472379.00,1756481.0,1756481.0] #-472379 vs. -473379
vT_1_2H=[14184.0,304.6,148.0,447.2,447.2]


EGG=0
LARVAE=1
PUPAE=2
ADULT1=3
ADULT2=4
WATER=5


#<precipitation related functionality v>


#NOTE:EPA evaporation in Technical guidance for hazards analysis, eq (7) page G-3
#TODO:check QS for spilled water
def QR(u,BS_s,T_t):#in l/day
    U=u * 1000.0/(3600.0)#step(wind_speed,t) * 1000.0/(60.0*60.0)#in m/s #km/h->m/s
    MW=18.01528#molecular wight of water in g/mol (gmol often called mol src:https://chemistry.stackexchange.com/questions/53455/g-gmol-vs-lb-lbmol)
    A=BS_s * 0.00107#in ft^2 #cm^2->ft^2
    VP=10.0**(8.07131-(1730.63/(233.426+T_t-273.15)))#Vapor pressure by Antoine Equation in mmHg
    R=82.05 #in atm cm 3 /g mole
    return ( (0.106 * U**0.78 * MW**(2.0/3.0)* A * VP)/(R* T_t) ) * 453.59237/(1.0/(60.0*24.0)) *1./1000. #Rate of release to air ml/day #(Ibs/min) ->  ml/day ->l/day

def QG(BS_s,p_t,t):#Quantity gathered#in litres
    return (BS_s * p_t*0.1) * 1.0 * 1./1000.#*1cm^3=1ml -> l

'''
    { QG(BS_s,t)-QR(BS_s,T(t))    if 0 < W < BS_c
dW= { QG(BS_s,t)                  if W <= 0.0
    { -QR(BS_s,T(t))               if W >= BS_c
Note: in the implementation we needed to add functions to make function continuous, otherwise odeint breaks
'''

def dW(W,BS_c,BS_s,T_t,p_t,wss_t,t):#in l/day
    epsilon=1e-3
    if(0+epsilon < W < BS_c-epsilon):
        return QG(BS_s,p_t,t)-QR(wss_t,BS_s,T_t)
    elif(W <= 0.0+epsilon):
        return QG(BS_s,p_t,t) - QR(wss_t,BS_s,T_t)*(W/epsilon)
    elif( W >= BS_c-epsilon):
        return QG(BS_s,p_t,t)*((BS_c-W)/epsilon) - QR(wss_t,BS_s,T_t)

def a0(W):
    return 70.0* W

def vGamma(vL,vBS_a,vW):
    return [gamma(vL[i],vBS_a[i],vW[i]) for i in range(0,len(vL))]

#</precipitation related functionality v>


def R_D(stage,T_t):#day^-1
    R=1.987 # universal gas constant
    R_D_298K=vR_D_298K[stage]
    #ro_25=ro_25_[stage]
    deltaH_A=vDeltaH_A[stage]
    deltaH_H=vDeltaH_H[stage]
    T_1_2H=vT_1_2H[stage] # K
    return R_D_298K * (T_t/298.0) * math.exp( (deltaH_A/R)* ((1.0/298.0)- (1.0/T_t)) ) / ( 1.0+ math.exp( (deltaH_H/R)* ((1.0/T_1_2H)-(1.0/T_t)) ) )


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

def beta(vW,vBS_od,vBS_id):#~1 if we have water,~0 if we dont
    return np.dot(vBS_od,vW/(vW+1e-4))+ np.dot(vBS_id,np.ones(len(vBS_id)))#TODO:check!

def dE(E,L,A1,A2,vW,T_t,BS_a,vBS_ic,vBS_od,vBS_id,n,m):
    egn=63.0*beta(vW,vBS_od,vBS_id)#The amount of eggs goes to zero when vW goes to zero.In the 1-dimensional case. when W=1e-3, egn(1e-3)~egn(0.5)/2#TODO:instead of 1e-3, it whould be a value related to the min water needed to lay eggs
    me=0.01#mortality of the egg, for T in [278,303]
    elr=R_D(EGG,T_t)
    ovr1=R_D(ADULT1,T_t)
    ovr2=R_D(ADULT2,T_t)
    v1=np.ones((n))
    #  ( (1,1,..1) - vGamma(vL,vBS_a,vW) ) . vE = (1- γ(vL[1],vBS_a[1],vW[1],...) ) . vE= Σ (1- γ(vL[i],vBS_a[i],vW[i])) * vE[i]
    inh_o=np.dot(v1 -vGamma(L*vBS_od,BS_a*vBS_od,vW    ) , E*vBS_od )
    v1=np.ones((m))
    #  ( (1,1,..1) - vGamma(vL,vBS_a,vBS_ic) ) . vE = (1-γ(vL[1],vBS_a[1],vBS_ic[1]) ,... ) . vE= Σ (1- γ(vL[i],vBS_a[i],vBS_ic[i])) * vE[i]
    inh_i=np.dot(v1 -vGamma(L*vBS_id,BS_a*vBS_id,vBS_ic) , E*vBS_id )
    return egn*( ovr1 *A1  + ovr2* A2) - me *E - elr* (inh_o + inh_i )

def dL(E,L,vW,T_t,BS_a,vBS_ic,vBS_od,vBS_id,n,m):
    elr=R_D(EGG,T_t)
    lpr=R_D(LARVAE,T_t)
    ml=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#mortality of the larvae, for T in [278,303]
    alpha0=1.5#Parameter to be fitted #1.0#HARDCODED!!!
    alpha=alpha0/BS_a#Σ vα[i] * vL[i]^2= Σ α0/(BS_a* vBS_d[i]) * (L*vBS_d[i])^2 = Σ α0/BS_a * L^2 * vBS_d[i] = α *L^2 *Σ BS_d[i]=α L^2 #Note: on why this is ok even though BS changed.
    v1=np.ones((n))
    inh_o=np.dot(v1 -vGamma(L*vBS_od,BS_a*vBS_od,vW    ) , E*vBS_od )
    v1=np.ones((m))
    inh_i=np.dot(v1 -vGamma(L*vBS_id,BS_a*vBS_id,vBS_ic) , E*vBS_id )
    return elr* (inh_o+inh_i ) - ml*L - alpha* L*L - lpr *L -35.6464*(1.-beta(vW,vBS_od,vBS_id))*L#-24.*(1.-beta(vW))*L# -log(1e-4/5502.)/(1.)=17.823207313460703

def dP(L,P,T_t):
    lpr=R_D(LARVAE,T_t)
    par=R_D(PUPAE,T_t)
    mp=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#death of pupae
    return lpr*L - mp*P  - par*P

def dA1(P,A1,T_t):
    par=R_D(PUPAE,T_t)
    ovr1=R_D(ADULT1,T_t)
    ef=0.83#emergence factor
    ma=0.09#for T in [278,303]
    return par*ef*P/2.0 - ma*A1 - ovr1*A1

def dA2(A1,A2,T_t):
    ovr1=R_D(ADULT1,T_t)
    ma=0.09#for T in [278,303]
    return ovr1*A1 - ma*A2

def diff_eqs(Y,t,parameters):
    '''The main set of equations'''
    T_t=parameters.weather.T(t)
    p_t=parameters.weather.p(t)
    ws_t=parameters.weather.ws(t)
    wss_t=parameters.ws_s*ws_t
    BS_a,vBS_oc,vBS_ic,vBS_od,vBS_id,vBS_os,n,m=parameters.BS_a,parameters.vBS_oc,parameters.vBS_ic,parameters.vBS_od,parameters.vBS_id,parameters.vBS_os,parameters.n,parameters.m

    E,L,P,A1,A2=Y[:WATER]
    vW=np.array(Y[WATER:])

    dY=np.zeros((5+n))
    dY[EGG]    = dE(E,L,A1,A2,vW,T_t,BS_a,vBS_ic,vBS_od,vBS_id,n,m)
    dY[LARVAE] = dL(E,L,vW,T_t,      BS_a,vBS_ic,vBS_od,vBS_id,n,m)
    dY[PUPAE]  = dP(L,P,T_t)
    dY[ADULT1] = dA1(P,A1,T_t)
    dY[ADULT2] = dA2(A1,A2,T_t)
    dY[WATER:] = [dW(vW[i],vBS_oc[i],vBS_os[i],T_t,p_t,wss_t,t) for i in range(0,n)]

    return dY   # For odeint
