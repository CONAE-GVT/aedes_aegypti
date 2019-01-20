import math
import numpy as np
from equations import vR_D,dvP,dA1,dA2


def ovsp(vW,vBS_d,vW_l,mBS_l):#OViposition Site Preference
    epsilon=1e-4
    vf=vW/(vW+epsilon) * vBS_d#check this method is not spontaneus generation of eggs.(like inventing adults.)
    if(np.all(vf)>epsilon): vf/=vf.sum()
    return np.where(mBS_l==np.floor(vW_l),1,0)*vf

def wetMask(vW_l,mBS_l):
    mask=np.where(mBS_l<vW_l,1,0)
    return np.where(mBS_l==np.floor(vW_l),vW_l%1,mask)


def dvE(mE,A1,A2,vW_t,vBS_d,elr,ovr1,ovr2,wet_mask,vW_l,mBS_l):
    egn=63.0
    me=0.01#mortality of the egg, for T in [278,303]
    ovsp_t=ovsp(vW_t,vBS_d,vW_l,mBS_l)
    return egn*( ovr1 *A1  + ovr2* A2)*ovsp_t - me * mE - elr * mE*wet_mask

def dvL(mE,vL,vW,T_t,BS_a,vBS_d,elr,lpr,vAlpha0,wet_mask):
    ml=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#mortality of the larvae, for T in [278,303]
    vAlpha=vAlpha0/(BS_a*vBS_d)
    return elr * np.sum(mE*wet_mask,axis=0) - ml*vL - vAlpha* vL*vL - lpr *vL#-35.6464*(1.-beta(vW,vBS_od,vBS_id))*L#-24.*(1.-beta(vW))*L# -log(1e-4/5502.)/(1.)=17.823207313460703


def diff_eqs(Y,t,parameters):
    '''The main set of equations'''
    T_t=parameters.weather.T(t)
    elr,lpr,par,ovr1,ovr2=vR_D(T_t)
    BS_a,vBS_h,vBS_d,vAlpha0,n,BS_l,mBS_l=parameters.BS_a,parameters.vBS_h,parameters.vBS_d,parameters.vAlpha0,parameters.n,parameters.BS_l,parameters.mBS_l
    EGG,LARVAE,PUPAE,ADULT1,ADULT2=parameters.EGG,parameters.LARVAE,parameters.PUPAE,parameters.ADULT1,parameters.ADULT2

    vW_t=parameters.vW(t)
    mE,vL,vP,A1,A2=Y[EGG].reshape((parameters.BS_l,n)),Y[LARVAE],Y[PUPAE],Y[ADULT1],Y[ADULT2]
    vW_l=vW_t/vBS_h * BS_l
    wet_mask=wetMask(vW_l,mBS_l)

    dY=np.empty( Y.shape )
    dY[EGG]    = dvE(mE,A1,A2,vW_t,vBS_d,elr,ovr1,ovr2,wet_mask,vW_l,mBS_l).reshape((1,BS_l*n))
    dY[LARVAE] = dvL(mE,vL,vW_t,T_t,BS_a,vBS_d,elr,lpr,vAlpha0,wet_mask)
    dY[PUPAE]  = dvP(vL,vP,T_t,lpr,par)
    dY[ADULT1] = dA1(vP,A1,T_t,par,ovr1)
    dY[ADULT2] = dA2(A1,A2,T_t,ovr1)

    return dY   # For odeint
