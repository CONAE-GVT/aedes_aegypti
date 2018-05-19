from otero_precipitation import Model
from equations import vR_D,vGamma
from config import Configuration
from scipy import optimize
import numpy as np
import pylab as pl
import datetime
import utils
import math


def f_error(x,f,domain,real_values):
    error=real_values - [ f(x,v) for v in domain]
    print(np.dot(error,error))
    return np.dot(error,error)#for leastsq return just error

class LabEquations:
    def __init__(self):
        pass

    def dvL(self,vE,vL,vW,T_t,      BS_a,vBS_d,elr,lpr,vAlpha0,a,b,c):
        ml=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#mortality of the larvae, for T in [278,303]
        #vAlpha0=a*np.exp(b*(vL-c))#this one goes to hell!
        vAlpha0=a/(0.05+np.exp(-c*(vL-b)))#a/(b+np.exp(-c*(vL-62))) <- candidate
        return elr* (1-vGamma(vL,BS_a*vBS_d,vW)) * vE - ml*vL - lpr *vL - vAlpha0*vL

    def __call__(self,Y,t,parameters):
        T_t=parameters.weather.T(t)
        p_t=parameters.weather.p(t)
        RH_t=parameters.weather.RH(t)
        elr,lpr,par,ovr1,ovr2=vR_D(T_t)
        BS_a,vBS_oc,vBS_ic,vBS_d,vBS_os,vAlpha0,n,m=parameters.BS_a,parameters.vBS_oc,parameters.vBS_ic,parameters.vBS_d,parameters.vBS_os,parameters.vAlpha0,parameters.n,parameters.m
        EGG,LARVAE,PUPAE,ADULT1,ADULT2,WATER=parameters.EGG,parameters.LARVAE,parameters.PUPAE,parameters.ADULT1,parameters.ADULT2,parameters.WATER

        vE,vL,vP,A1,A2,vW=Y[EGG],Y[LARVAE],Y[PUPAE],Y[ADULT1],Y[ADULT2],Y[WATER]
        vW=np.concatenate((vW,vBS_ic))#constant water for inside BS
        a,b,c=parameters.a,parameters.b,parameters.c
        dY=np.zeros((3*(n+m)+2+n))
        dY[EGG]    = np.array([0]*(n+m))#no eggs
        dY[LARVAE] = self.dvL(vE,vL,vW,T_t,      BS_a,vBS_d,elr,lpr,vAlpha0,a,b,c)
        dY[PUPAE]  = lpr*vL#no mortality and never becoming adult
        dY[ADULT1] = np.array([0]*(n+m))#no adults
        dY[ADULT2] = np.array([0]*(n+m))#no adults
        dY[WATER] = [self.dW(vW[i],vBS_oc[i],vBS_os[i],T_t,p_t,RH_t,t) for i in range(0,n)]#note that this goes till n.

        return dY

def macia(x,v):
    a,b,c=x[0],x[1],x[2]
    L0=v
    W,t_final=0.175,44.
    configuration=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'amount':1,
            'outside_capacity':[],
            'outside_surface':[],
            'outside_distribution':[],
            'inside_distribution':[1.],
            'inside_capacity':[W]
            },
        'simulation':{
            'start_date':datetime.date(2017,3,1),
            'end_date':datetime.date(2017,4,15),
            'initial_condition':[0.]*1 + [L0]*1 +[0.]*1 + [0.,0.]
            },
        'biology':{
            'alpha0':[0]#alpha is not used
            }
        })

    model=Model(configuration)
    #lab environment
    T=lambda t: 26.+273.15#26 +/- 2
    model.parameters.weather.T=T
    model.parameters.a,model.parameters.b,model.parameters.c=a,b,c
    time_range,INPUT,RES=model.solveEquations(equations=LabEquations())
    LARVAE,PUPAE,WATER=model.parameters.LARVAE,model.parameters.PUPAE,model.parameters.WATER
    vBS_ic=model.parameters.vBS_ic
    L_initial,P_final=RES[0,LARVAE],RES[-1,PUPAE]
    return (P_final/L_initial)[0]


def fit(x0,f,bounds,domain,real_values,method='SLSQP'):
    #Find optimal parameters
    res=None
    if(method=='SLSQP'):
        res=optimize.minimize(f_error,x0,(f,domain,real_values),method='SLSQP',bounds=bounds)#,constraints=constraints
    else:
        res=optimize.differential_evolution(f_error,bounds=bounds,args=(domain,real_values))

    print(res)
    x=res.x
    pl.plot(range(0,domain.max()),[f(x,v) for v in range(0,domain.max())],label=f.__name__ )
    pl.plot(domain,real_values, '^y',label='')
    pl.xlabel('')
    pl.ylabel('')
    pl.legend(loc=0)

    pl.show()

def fitMacia():
    #define values to fit
    domain=     np.array([4,8   ,16  ,32  ,64    ,128  ,256])
    survival=np.array([1,1   ,1   ,1   ,0.970 ,0.450,0.065])
    #find optimal
    bounds=((0,10),(54,74),(0,10))#((0,5),(0.05,5),(0,5))<- candidate
    x0=np.array([0,0,0])
    fit(x0,macia,bounds,domain,survival)

if(__name__ == '__main__'):
    fitMacia()
