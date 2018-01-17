#coding: utf-8
from scipy import interpolate
import scipy.integrate as spi
import numpy as np
import pylab as pl
import datetime
import fourier
import fitter
import utils
import math
import rk



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

BS_a=50.0#amount of BS #1.0/21.0#??
BS_o=0.5#1.0/2.0#proportion of BS that are affeted by rain (are outside)
vBS_oc, vBS_ic=np.array([.04,0.5,20.]), np.array([0.5])#[1./2.]#capacity in litres
vBS_od, vBS_id=BS_o*np.array([0.1,0.4,0.5]), (1.-BS_o)*np.array([1.])#distribution of BS, the sum must be equal to 1.0#[1.0]
vBS_os, vBS_is=math.pi*np.array([1.25**2,5.25**2,15.25**2]), math.pi*np.array([5.25**2])#in cm^2 #TODO:CHECK THAT ALL UNITS ARE CONSISTENT!!!
n,m=len(vBS_oc),len(vBS_ic)
ws_s=0.5#wind shield in [0,1]
#a0=70.0* BS_capacity #??
#Cordoba
location,start_date,end_date={'name':'cordoba','station':'SACO','zones':['Centro','NO','NE','SE','SO']},datetime.date(2014, 7, 1),datetime.date(2017, 10, 1)
#Jujuy
#location,start_date,end_date={'name':'jujuy libertador','station':'SASJ','zones':[]},datetime.date(2010, 07, 1),datetime.date(2013, 7, 1)


AEDIC_INDICES_FILENAME='data/private/Indices aedicos Historicos '+location['name']+'.xlsx'
WEATHER_STATION_DATA_FILENAME='data/public/wunderground_'+location['station']+'.csv'

#<precipitation related functionality v>
'''
Convenience method
Let  g=getAsLambdaFunction(aps2,precipitations)
calling g(t) would be as calling aps2(precipitations,t)
'''
def getAsLambdaFunction(f,values):
    return lambda t: f(values,t)

#Area preserving Sin
def aps(values,t):
    return values[int(t)] * (math.sin(2.*math.pi*t + 3.*math.pi/2.) +1.)

#Step
def step(values,t):
    return values[int(t)]

precipitations = utils.getPrecipitationsFromCsv(WEATHER_STATION_DATA_FILENAME,start_date,end_date)
p=getAsLambdaFunction(aps,precipitations)

wind_speed=utils.getMeanWindSpeedFromCsv(WEATHER_STATION_DATA_FILENAME,start_date,end_date)
ws=fourier.fourier(wind_speed,50)

#NOTE:EPA evaporation in Technical guidance for hazards analysis, eq (7) page G-3
#TODO:check QS for spilled water
def QR(u,BS_s,T_t):#in l/day
    U=u * 1000.0/(3600.0)#step(wind_speed,t) * 1000.0/(60.0*60.0)#in m/s #km/h->m/s
    MW=18.01528#molecular wight of water in g/mol (gmol often called mol src:https://chemistry.stackexchange.com/questions/53455/g-gmol-vs-lb-lbmol)
    A=BS_s * 0.00107#in ft^2 #cm^2->ft^2
    VP=10.0**(8.07131-(1730.63/(233.426+T_t-273.15)))#Vapor pressure by Antoine Equation in mmHg
    R=82.05 #in atm cm 3 /g mole
    return ( (0.106 * U**0.78 * MW**(2.0/3.0)* A * VP)/(R* T_t) ) * 453.59237/(1.0/(60.0*24.0)) *1./1000. #Rate of release to air ml/day #(Ibs/min) ->  ml/day ->l/day

def QG(BS_s,t):#Quantity gathered#in litres
    return (BS_s * p(t)*0.1) * 1.0 * 1./1000.#*1cm^3=1ml -> l

'''
    { QG(BS_s,t)-QR(BS_s,T(t))    if 0 < W < BS_c
dW= { QG(BS_s,t)                  if W <= 0.0
    { -QR(BS_s,T(t))               if W >= BS_c
Note: in the implementation we needed to add functions to make function continuous, otherwise odeint breaks
'''

def dW(W,BS_c,BS_s,T_t,t):#in l/day
    epsilon=1e-3
    if(0+epsilon < W < BS_c-epsilon):
        return QG(BS_s,t)-QR(ws_s*ws(t),BS_s,T_t)
    elif(W <= 0.0+epsilon):
        return QG(BS_s,t) - QR(ws_s*ws(t),BS_s,T_t)*(W/epsilon)
    elif( W >= BS_c-epsilon):
        return QG(BS_s,t)*((BS_c-W)/epsilon) - QR(ws_s*ws(t),BS_s,T_t)

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


#a,b,c = fitter.getOptimalParameters(utils.getAverageTemperaturesFromCsv(WEATHER_STATION_DATA_FILENAME,start_date,end_date))#parameters for T(t)
#def T(t):#K, beginning on the first of july
#    return a + b * math.cos( (2.0* math.pi * t)/365.25 + c) + 273.15 #To Kelvin
#T=fourier.fourier(utils.getAverageTemperaturesFromCsv(WEATHER_STATION_DATA_FILENAME,start_date,end_date),50)
#time_domain,min_max_temperatures=utils.getMinMaxTemperaturesFromCsv(WEATHER_STATION_DATA_FILENAME,start_date,end_date)
#T = interpolate.interp1d(time_domain,min_max_temperatures,'cubic')
T=interpolate.InterpolatedUnivariateSpline(range(0,(end_date - start_date).days),utils.getAverageTemperaturesFromCsv(WEATHER_STATION_DATA_FILENAME,start_date,end_date))

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

def beta(vW):#~1 if we have water,~0 if we dont
    return np.dot(vBS_od,vW/(vW+1e-4))+ np.dot(vBS_id,np.ones(len(vBS_id)))#TODO:check!

def dE(E,L,P,A1,A2,vW,T_t):
    egn=63.0*beta(vW)#The amount of eggs goes to zero when vW goes to zero.In the 1-dimensional case. when W=1e-3, egn(1e-3)~egn(0.5)/2#TODO:instead of 1e-3, it whould be a value related to the min water needed to lay eggs
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

def dL(E,L,P,A1,A2,vW,T_t):
    elr=R_D(EGG,T_t)
    lpr=R_D(LARVAE,T_t)
    ml=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#mortality of the larvae, for T in [278,303]
    alpha0=1.5#Parameter to be fitted #1.0#HARDCODED!!!
    alpha=alpha0/BS_a#Σ vα[i] * vL[i]^2= Σ α0/(BS_a* vBS_d[i]) * (L*vBS_d[i])^2 = Σ α0/BS_a * L^2 * vBS_d[i] = α *L^2 *Σ BS_d[i]=α L^2 #Note: on why this is ok even though BS changed.
    v1=np.ones((n))
    inh_o=np.dot(v1 -vGamma(L*vBS_od,BS_a*vBS_od,vW    ) , E*vBS_od )
    v1=np.ones((m))
    inh_i=np.dot(v1 -vGamma(L*vBS_id,BS_a*vBS_id,vBS_ic) , E*vBS_id )
    return elr* (inh_o+inh_i ) - ml*L - alpha* L*L - lpr *L -35.6464*(1.-beta(vW))*L#-24.*(1.-beta(vW))*L# -log(1e-4/5502.)/(1.)=17.823207313460703

def dP(E,L,P,A1,A2,T_t):
    lpr=R_D(LARVAE,T_t)
    par=R_D(PUPAE,T_t)
    mp=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#death of pupae
    return lpr*L - mp*P  - par*P

def dA1(E,L,P,A1,A2,T_t):
    par=R_D(PUPAE,T_t)
    ovr1=R_D(ADULT1,T_t)
    ef=0.83#emergence factor
    ma=0.09#for T in [278,303]
    return par*ef*P/2.0 - ma*A1 - ovr1*A1

def dA2(E,L,P,A1,A2,T_t):
    ovr1=R_D(ADULT1,T_t)
    ma=0.09#for T in [278,303]
    return ovr1*A1 - ma*A2

def diff_eqs(V,t):
    '''The main set of equations'''
    Y=np.zeros((5+n))
    #V[V<0]=0#this is to make rk work
    T_t=T(t)
    vW=np.array([V[WATER+i] for i in range(0,n)])
    Y[EGG] = dE(V[0],V[1],V[2],V[3],V[4],vW,T_t)
    Y[LARVAE] = dL(V[0],V[1],V[2],V[3],V[4],vW,T_t)
    Y[PUPAE] = dP(V[0],V[1],V[2],V[3],V[4],T_t)
    Y[ADULT1] = dA1(V[0],V[1],V[2],V[3],V[4],T_t)
    Y[ADULT2] = dA2(V[0],V[1],V[2],V[3],V[4],T_t)
    for i in range(0,n):
        Y[WATER+i]=dW(V[5+i],vBS_oc[i],vBS_os[i],T_t,t)

    return Y   # For odeint

def solveEquations(INPUT = [100.0, 0.0,0.0,0.0,0.0]+ [0. for i in range(0,n)]):
    time_range = np.linspace(0, (end_date - start_date).days-2, (end_date - start_date).days * 10)
    RES = spi.odeint(diff_eqs,INPUT,time_range,hmax=0.01)
    #RES=rk.solve(diff_eqs,INPUT,time_range)
    #RES=rk.scipy_solve(diff_eqs,INPUT,time_range,'dopri',{'max_step':time_range[1]-time_range[0],'rtol':1e-3, 'atol':1e-6} )
    #RES = spi.odeint(op.diff_eqs,INPUT,time_range,hmax=0.5,rtol=[1e-2]*5 +[1e-2]*op.n ,atol=[1]*5 +[1e-4]*op.n)#,hmax=0.01
    return time_range,INPUT,RES


if(__name__ == '__main__'):
    #Solve differential equations
    time_range,INPUT,RES=solveEquations()

    #Ploting
    #Amount of larvaes,pupaes and adults
    pl.subplot(411)
    pl.plot(time_range,RES[:,LARVAE], '-r', label='L')
    pl.plot(time_range,RES[:,PUPAE], '-g', label='P')
    pl.plot(time_range,RES[:,ADULT1], '-b', label='A1')
    pl.plot(time_range,RES[:,ADULT2], '-m', label='A2')
    pl.xlabel('Time(in days starting in July)')
    pl.ylabel('')
    pl.legend(loc=0)

    #Temperature
    vT=np.vectorize(T)
    pl.subplot(412)
    #pl.plot(time_range,vT(time_range), '-k', label='Temperature')
    for i in range(0,n):
        pl.plot(time_range,RES[:,WATER+i], label='W for %sL, %scm^2'%(vBS_oc[i],vBS_os[i]) )
    pl.xlabel('Time(in days starting in July)')
    #pl.ylabel('Temperature(in K)')
    pl.legend(loc=0)

    #Indices
    pl.subplot(414)
    hi_days, his=utils.getIndexesForPlot(AEDIC_INDICES_FILENAME,start_date,0,3)
    pl.plot(hi_days,np.multiply(1.0/max(his),his), '^m', label='House Indices normalized',clip_on=False, zorder=100,markersize=8) #BEWARE!!! they are all relative, so no comparison in absolute terms
    bi_days, bis=utils.getIndexesForPlot(AEDIC_INDICES_FILENAME,start_date,0,4)
    pl.plot(bi_days,np.multiply(1.0/max(bis),bis), '^c', label='Breteau Indices normalized',clip_on=False, zorder=100,markersize=8)
    ri_days, ris=utils.getIndexesForPlot(AEDIC_INDICES_FILENAME,start_date,0,5)
    pl.plot(ri_days,np.multiply(1.0/max(ris),ris), '^y', label='Recip. Indices normalized',clip_on=False, zorder=100,markersize=8)
    pl.plot(time_range, np.multiply(1.0/max(RES[:,LARVAE]), RES[:,LARVAE]), '-y', label='Larvaes normalized')
    pl.xlabel('Time(in days starting in July)')
    locations,labels= utils.getDatesTicksForPlot(start_date, time_range)
    pl.xticks(locations,labels,rotation='vertical')
    pl.ylabel('')
    pl.grid()
    pl.legend(loc=0,prop={'size':10})

    #indices by zone
    pl.subplot(413)
    for zone in location['zones']:
        hi_days, his= utils.getIndexesForPlot(AEDIC_INDICES_FILENAME,start_date,0,2,zone,1)
        pl.plot(hi_days,np.multiply(1.0/max(his),his), '^', label=zone +' normalized',clip_on=False, zorder=100,markersize=8) #BEWARE!!! they are all relative, so no comparison in absolute terms
    pl.xlabel('Time(in days starting in July)')
    pl.ylabel('')
    pl.legend(loc=0,prop={'size':10})

    pl.show()
