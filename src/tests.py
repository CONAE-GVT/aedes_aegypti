import math
import utils
import datetime
import pylab as pl
import numpy as np
import matplotlib
import matplotlib.dates
from scipy.stats import stats
import otero_precipitation as op

def discreteQG(BS_s,precipitations,t):#l/day
    return (BS_s * precipitations[int(t)]*0.1) * 1.0 * 1./1000.#*1cm^3=1ml -> litres

def discreteW(a0,BS_c,BS_s,precipitations):#in l
    a=[0]* len(precipitations)
    a[0]=a0
    #here we are counting the amount rained, at the END of the day. That's why we need use the t-1 instead of just t
    for t in range(1,len(precipitations)):
        a[t]=max(min(a[t-1] + discreteQG(BS_s,precipitations,t-1) - op.QR(op.ws_s*op.ws(t-1),BS_s,op.T(t-1)),BS_c),0.0)
    return a

def step(values,t):
    return values[int(t)]

def printCorrelation():
    time_range,INPUT,RES=op.solveEquations()
    vT=np.vectorize(op.T)
    print('(Pearson s correlation coefficient,2-tailed p-value): ' + str(stats.pearsonr(RES[:,op.LARVAE], vT(time_range))) )

#tests
def testAps():
    days=range(0,(op.end_date - op.start_date).days)
    s=op.getAsLambdaFunction(step,op.precipitations)
    print('aps_test')
    for t in days:
        aps_test=spi.quad(op.p,t,t+1)
        step_test=spi.quad(s,t,t+1)
        if(abs(op.precipitations[t]-aps_test[0])>1e-10):
            print("%i: %f vs. %f " % (t,aps_test[0],op.precipitations[t]) )
            #print("%f: %f vs. %f " % (t,aps_test[0],step_test[0]))
def testdW():
    time_range,INPUT,RES=op.solveEquations()
    days=range(0,(op.end_date - op.start_date).days)
    accumulated_water=[discreteW(INPUT[op.WATER+i],op.vBS_oc[i],op.vBS_os[i],op.precipitations) for i in range(0,op.n)]
    print('dW_test')
    #compare
    for i in range(0,op.n):
        for t in days:
            if(abs(accumulated_water[i][t]-RES[t,op.WATER+i])>5e-2):
                print("%i: %f - %f = %f" % (t,RES[t,op.WATER+i],accumulated_water[i][t],RES[t,op.WATER+i]-accumulated_water[i][t]) )
    #plot
    for i in range(0,op.n):
        pl.plot(range(0,len(accumulated_water[i])),accumulated_water[i], label='W (Discrete version) for %sL, %scm^2'%(op.vBS_oc[i],op.vBS_os[i]) )
        pl.plot(time_range,RES[:,op.WATER+i], label='W (Continuous version) for %sL, %scm^2'%(op.vBS_oc[i],op.vBS_os[i]) )
    #pl.plot(range(0,600),[gamma(L,0.1) for L in range(0,600)], '-r', label='gamma')
    pl.xlabel('Time(in days starting in July)')
    pl.ylabel('Water in Litres')
    pl.legend(loc=0)
    pl.show()

def testRESvsOldRES(filename):
    time_range,INPUT,RES=op.solveEquations()
    print("||old_RES-RES|| = "),
    old_RES=np.load(filename)
    print(np.linalg.norm(old_RES-RES))
    #np.save('backup/RES_%s.npy'%(datetime.datetime.now().isoformat()),RES)


def testModel(vBS_oc=np.array([0.5]),vBS_od=np.array([1.]), precipitations=op.precipitations):
    assert(np.sum(vBS_od)==1)
    op.BS_o=1.#proportion of BS that are affected by rain (are outside)
    op.vBS_oc=vBS_oc#[1./2.]#capacity in litres
    op.vBS_od,op.vBS_id=op.BS_o*vBS_od, (1.-op.BS_o)*np.array([1.])#distribution of BS, the sum must be equal to 1.0#[1.0]
    op.vBS_os=math.pi*np.array([5.25**2]*len(op.vBS_oc))#in cm^2
    op.n=len(op.vBS_oc)

    op.precipitations=precipitations
    op.p=op.getAsLambdaFunction(op.aps,op.precipitations)
    time_range,INPUT,RES=op.solveEquations([100.0, 0.0,0.0,0.0,0.0]+ [0./2. for i in range(0,op.n)])

    #Ploting
    pl.figure()
    #Amount of larvaes,pupaes and adults
    ax1 =  pl.subplot(311)
    date_range=[datetime.timedelta(days=d)+datetime.datetime.combine(op.start_date,datetime.time()) for d in time_range]
    ax1.xaxis.set_major_locator( matplotlib.dates.MonthLocator() )
    ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%Y-%m-%d') )
    ax1.xaxis.set_minor_locator( matplotlib.dates.DayLocator() )

    pl.plot(date_range,RES[:,op.EGG]*5e-2, '-k', label='E  5e-2')
    pl.plot(date_range,RES[:,op.LARVAE], '-r', label='L')
    pl.plot(date_range,RES[:,op.PUPAE], '-g', label='P')
    pl.plot(date_range,RES[:,op.ADULT1], '-b', label='A1')
    pl.plot(date_range,RES[:,op.ADULT2], '-m', label='A2')
    #pl.plot(date_range,RES[:,op.ADULT2]+RES[:,op.ADULT1], '-m', label='A1+A2')
    pl.xlabel('Time(in days starting in July)')
    pl.ylabel('')
    pl.legend(loc=0)
    pl.xticks(rotation='vertical')

    #Water in containers(in L)
    pl.subplot(312,sharex=ax1)#sharex to zoom all subplots if one is zoomed
    for i in range(0,op.n):
        pl.plot(date_range,RES[:,op.WATER+i], label='W(t) for %sL, %scm^2'%(op.vBS_oc[i],op.vBS_os[i]) )
    pl.xlabel('Time(in days starting in July)')
    pl.ylabel('Litres')
    pl.legend(loc=0)
    pl.xticks(rotation='vertical')

    #precipitations(in mm.)
    pl.subplot(313,sharex=ax1)
    pl.plot(date_range,[op.p(t) for t in time_range],'-b', label='p(t)')
    pl.xlabel('Time(in days starting in July)')
    pl.ylabel('mm.')
    pl.legend(loc=0)
    pl.xticks(rotation='vertical')

    getLarvaeSurvivalDry(time_range,RES[:,op.WATER:op.WATER+op.n],RES[:,op.LARVAE],{'start':100,'end':1600})
    #plotBeta()

def getLarvaeSurvivalDry(time_range,vW,L,time_window):
    #print(W[:,0]<1e-3)
    #print(time_range[W<1e-3])
    epsilon=1e-4
    times_dry_container=[t for i,t in enumerate(time_range) if np.dot(op.vBS_od,vW[i,:])<epsilon and time_window['start']<t<time_window['end'] ]
    times_death_larvae= [t for i,t in enumerate(time_range) if L[i]<epsilon and time_window['start']<t<time_window['end']]
    if not times_dry_container or not times_death_larvae:
        return
    print('W<epsilon at t=%f with epsilon=%f'%(times_dry_container[0],epsilon) )
    print('L<epsilon at t=%f with epsilon=%f'%(times_death_larvae[0],epsilon) )
    print(L[np.where(time_range==times_dry_container[0])])

def plotBeta():
    pl.figure()
    W_range=np.linspace(5, 0,5*20)
    pl.plot(W_range,[op.beta(np.array([W])) for W in W_range])

#TODO:check this method!!!
def getRainChances(d):
    rainy_months=[1,2,3,10,11,12]
    if (datetime.timedelta(days=d)+op.start_date).month in rainy_months:
        return 0.3
    else:
        return 0.005
def getFakePrecipitation(days):
    total_precipitation=1000.#mm
    tmp=np.array([np.random.rand(1) if np.random.rand(1)<getRainChances(d) else 0. for d in range(0,days)])#np.random.normal(50,15)
    return (total_precipitation/np.sum(tmp))*tmp

if(__name__ == '__main__'):
    #printCorrelation()
    #testAps()
    #testdW()
    testRESvsOldRES('odeint_otero_precipitation.npy')
    quit()

    vvBS_oc=np.array([np.array([0.1]),np.array([0.5]),np.array([1.2])])#in litres
    #pathological fake precipitations

    #for vBS_oc in vvBS_oc:
    #    testModel(vBS_oc=vBS_oc,precipitations=[0 if d>500 else 15. for d in range(0,(op.end_date - op.start_date).days)])
    #    testModel(vBS_oc=vBS_oc,precipitations=[0. if 365<d<600 else 15. for d in range(0,(op.end_date - op.start_date).days)])
    #    testModel(vBS_oc=vBS_oc,precipitations=[6. for d in range(0,(op.end_date - op.start_date).days)])#I should remove the wind.Maybe using Ivanov model (Romanenko 1961)

    #fake precipitations
    fake_precipitations=getFakePrecipitation((op.end_date - op.start_date).days)
    #for vBS_oc in vvBS_oc:
    #    testModel(vBS_oc=vBS_oc,precipitations=fake_precipitations)
    #testModel(vBS_oc=vvBS_oc.flatten(),vBS_od=np.array([0.1,0.3,0.6]),precipitations=fake_precipitations)
    #real precipitations
    #for vBS_oc in vvBS_oc:
    #    testModel(vBS_oc=vBS_oc)
    testModel(vBS_oc=np.array([0.5]),vBS_od=np.array([1.]))
    pl.show()
