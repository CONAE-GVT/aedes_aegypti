from equations import *

def dvE(vE,F,vBS_d,ovr,elr):
    egn=63.0
    me=0.01#mortality of the egg, for T in [278,303]
    return egn* np.expand_dims(ovr * F,axis=2)*vBS_d - (me+elr) * vE

def dvL(vE,vL,T_t,elr,lpr,vAlpha):
    ml=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#mortality of the larvae, for T in [278,303]
    return elr * vE - (ml+ vAlpha* vL+ lpr) *vL

def dvP(vL,vP,T_t,lpr,par):
    mp=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#death of pupae
    return lpr*vL - (mp+par)*vP

def dA1(vP,A1,par,cycle1):
    ef=0.83#emergence factor
    ma=0.091#for T in [278,303]
    return np.sum(par*ef*vP/2.0,axis=2) - ma*A1 - cycle1*A1

def dF(A1,F,A2,P,ovr,cycle1,cycle2):
    ma=0.091#for T in [278,303]
    beta_p=830./100.**2#TODO: check if its correct!!!# dispersal coefficient for perpendicular flights
    beta_d=beta_p/2.#dispersal coefficient for diagonal flights
    return cycle1*A1 + cycle2*A2 - (ovr + ma + 4.*beta_d + 4*beta_p)*F +\
        beta_p*(np.roll(F,(0, 1),axis=(1,0))*P[:,:,0] +  np.roll(F,(0,-1),axis=(0,1))*P[:,:,2] +  np.roll(F,(-1,0),axis=(0,1))*P[:,:,4] + np.roll(F,(1,0),axis=(1,0))*P[:,:,6] ) +\
        beta_d*(np.roll(F,(-1, 1),axis=(1,0))*P[:,:,1] + np.roll(F,(-1,-1),axis=(1,0))*P[:,:,3] + np.roll(F,(1,-1),axis=(1,0))*P[:,:,5] + np.roll(F,(1,1),axis=(1,0))*P[:,:,7] )#TODO:Check if this is correct!

def dA2(F,A2,ovr,cycle2):
    ma=0.091#for T in [278,303]
    return ovr*F - cycle2*A2 - ma*A2

def diff_eqs(Y,t,parameters):
    '''The main set of equations'''
    T_t=parameters.weather.T(t)
    elr,lpr,par,cycle1,cycle2=vR_D(T_t)
    vBS_a,vBS_d,vAlpha,ovr,P,n=parameters.vBS_a,parameters.vBS_d,parameters.vAlpha,parameters.ovr,parameters.P,parameters.n
    HEIGHT,WIDTH=P.shape[:2]
    EGG,LARVAE,PUPAE,ADULT1,FLYER,ADULT2=parameters.EGG,parameters.LARVAE,parameters.PUPAE,parameters.ADULT1,parameters.FLYER,parameters.ADULT2

    Y=Y.reshape(HEIGHT,WIDTH,3*n + 3)
    vW_t=parameters.vW(t)
    vE,vL,vP,A1,F,A2=Y[:,:,EGG],Y[:,:,LARVAE],Y[:,:,PUPAE],Y[:,:,ADULT1],Y[:,:,FLYER],Y[:,:,ADULT2]

    dY=np.empty(Y.shape)
    dY[:,:,EGG]    = dvE(vE,F,vBS_d,ovr,elr)
    dY[:,:,LARVAE] = dvL(vE,vL,T_t,elr,lpr,vAlpha )
    dY[:,:,PUPAE]  = dvP(vL,vP,T_t,lpr,par)
    dY[:,:,ADULT1] = dA1(vP,A1,par,cycle1)
    dY[:,:,FLYER]  = dF(A1,F,A2,P,ovr,cycle1,cycle2)
    dY[:,:,ADULT2] = dA2(F,A2,ovr,cycle2)

    return dY.reshape(HEIGHT*WIDTH*(3*n + 3) )   # For odeint
