import math
import utils
from utils import getSurface,getCapacity#not sure if this is a good practice
import datetime
import pylab as pl
import numpy as np
import matplotlib
import matplotlib.dates
from scipy.stats import stats
from scipy import interpolate
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

a_time_range = np.linspace(0, (op.end_date - op.start_date).days-2, (op.end_date - op.start_date).days * 10)#TODO:find a better way...
calls=np.array([0]*len(a_time_range))
negatives=np.array([0]*len(a_time_range))
def decoratedEquations(Y,t):
    dY=op.diff_eqs(Y,t)
    calls[(np.abs(a_time_range-t)).argmin()]+=1
    if(np.any(Y<0)): negatives[(np.abs(a_time_range-t)).argmin()]+=1
    return dY

def normalize(values):#TODO:change name.
    return (values-np.min(values))/(np.max(values)-np.min(values))

def subData(time_range,Y,date_range,an_start_date):
    #get the index of an_start_date
    index=None
    for i,date in enumerate(date_range):
        if(date.date()==an_start_date):
            index=i
            break#conserve the first one.
    return time_range[index:],Y[index:,:],date_range[index:]

def testModel(BS_o=1.0,vBS_oc=np.array([0.5]),vBS_os=math.pi*np.array([5.25**2]),vBS_od=np.array([1.]), precipitations=op.precipitations,subplots=[['E','L'],['W']],plot_start_date=None):
    assert(np.sum(vBS_od)==1)
    op.BS_o=BS_o
    op.vBS_oc=vBS_oc#[1./2.]#capacity in litres
    op.vBS_od,op.vBS_id=op.BS_o*vBS_od, (1.-op.BS_o)*np.array([1.])#distribution of BS, the sum must be equal to 1.0#[1.0]
    op.vBS_os=vBS_os#in cm^2
    op.n=len(op.vBS_oc)

    op.precipitations=precipitations
    op.p=op.getAsLambdaFunction(op.aps,op.precipitations)
    time_range,INPUT,RES=op.solveEquations(initial_condition=[100.0, 0.0,0.0,0.0,0.0]+ [0./2. for i in range(0,op.n)],equations=decoratedEquations)

    #Ploting
    pl.figure()
    pl.subplots_adjust(top=0.95,hspace=0.28)
    ax1=None
    for i,subplot in enumerate(subplots):
        subplot_id=len(subplots)*100 + 10 + (i+1)
        if(i==0):
            ax1=pl.subplot(subplot_id)
            date_range=[datetime.timedelta(days=d)+datetime.datetime.combine(op.start_date,datetime.time()) for d in time_range]
            ax1.xaxis.set_major_locator( matplotlib.dates.MonthLocator() )
            ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%Y-%m-%d') )
            ax1.xaxis.set_minor_locator( matplotlib.dates.DayLocator() )
        else:
            pl.subplot(subplot_id,sharex=ax1)#sharex to zoom all subplots if one is zoomed

        if(plot_start_date):
            time_range,RES,date_range=subData(time_range,RES,date_range,plot_start_date)

        #Amount of larvaes,pupaes and adults
        if ('E' in subplot): pl.plot(date_range,RES[:,op.EGG]*5e-2, '-k', label='E *  5e-2')
        if ('L' in subplot): pl.plot(date_range,RES[:,op.LARVAE], '-r', label='L')
        if ('P' in subplot): pl.plot(date_range,RES[:,op.PUPAE], '-g', label='P')
        if ('A1' in subplot): pl.plot(date_range,RES[:,op.ADULT1], '-b', label='A1')
        if ('A2' in subplot): pl.plot(date_range,RES[:,op.ADULT2], '-m', label='A2')
        if ('A1+A2' in subplot): pl.plot(date_range,RES[:,op.ADULT2]+RES[:,op.ADULT1], '-m', label='A1+A2')
        if ('nLvsI' in subplot):
            ri_days, ris=utils.getIndexesForPlot(op.AEDIC_INDICES_FILENAME,op.start_date,0,5)
            pl.plot([datetime.timedelta(days=d)+datetime.datetime.combine(op.start_date,datetime.time()) for d in ri_days], np.array(ris)/max(ris), '^y', label='Recip. Indices normalized',clip_on=False, zorder=100,markersize=8)
            pl.plot(date_range,RES[:,op.LARVAE]/max(RES[:,op.LARVAE]), '-r', label='Larvaes normalized')
        if('nEvsO' in subplot):
            for i in range(1,30):
                ovitrap_eggs=utils.getOvitrapEggsFromCsv('data/private/Datos sensores de oviposicion.NO.csv',op.start_date,op.end_date,i)
                pl.plot([datetime.timedelta(days=d)+datetime.datetime.combine(op.start_date,datetime.time()) for d in range(0,len(ovitrap_eggs))], [e/max(ovitrap_eggs) if e else None for e in ovitrap_eggs], '^', label='Ovitrap %s eggs'%i,clip_on=False, zorder=100,markersize=8)
            pl.plot(date_range,RES[:,op.EGG]/max(RES[:,op.EGG]), '-k', label='Eggs normalized')
        if('lwE' in subplot):
            pl.plot(date_range, [RES[(np.abs(time_range-t)).argmin(),op.EGG]-RES[(np.abs(time_range-(t-7))).argmin(),op.EGG] for t in time_range], '-', label='Eggs oviposited this week')
        pl.ylabel('')

        #Complete lifecycle
        if('clc' in subplot):
            pl.plot(date_range,[sum([1./op.R_D(stage,op.T(t)) for stage in [op.EGG,op.LARVAE,op.PUPAE,op.ADULT1,op.ADULT2]]) for  t in time_range],label='Complete life cicle(from being an egg to the second oviposition)')
            pl.ylabel('days')
        #Water in containers(in L)
        if ('W' in subplot):
            for i in range(0,op.n):
                pl.plot(date_range,RES[:,op.WATER+i], label='W(t) for %sL, %scm^2, %s%%'%(op.vBS_oc[i],op.vBS_os[i],op.vBS_od[i]*100.) )
            pl.ylabel('Litres')

        #spaa vs cimsim
        if ('spaavscimsim' in subplot):
            for i in range(0,op.n):
                pl.plot(time_range,RES[:,op.WATER+i]*1000.0/op.vBS_os[i], label='W(t) for %sL, %scm^2, %s%%'%(op.vBS_oc[i],op.vBS_os[i],op.vBS_od[i]*100.) )#L->ml->mm->cm
            pl.plot(utils.getValuesFromCsv('data/test/cimsim_containers_2015_se09.csv',op.start_date,op.end_date,1,verbose=False),label='CIMSiM')

        #Temperature in K
        if ('T' in subplot):
            pl.plot(date_range,[op.T(t) for t in time_range], label='Temperature')
            pl.ylabel('K')

        #precipitations(in mm.)
        if ('p' in subplot):
            pl.plot(date_range,[op.p(t) for t in time_range],'-b', label='p(t)')
            pl.ylabel('mm./day')

        #Wind Speed(in km/h.)
        if ('ws' in subplot):
            pl.plot(date_range,[op.ws(t) for t in time_range], label='ws(t)')
            pl.ylabel('km/h')

        #Beta
        if ('b' in subplot):
            pl.plot(date_range,[op.beta(RES[(np.abs(time_range-t)).argmin(),op.WATER:]) for t in time_range], label='beta(vW)')
            pl.ylabel('')

        #debugging plots
        #Calls
        if ('c' in subplot):
            pl.plot(date_range,calls, label='calls')
            pl.ylabel('')
            print('Calls: %s'%sum(calls))
        #Negatives
        if ('n' in subplot):
            pl.plot(date_range,negatives, label='negatives')
            pl.ylabel('')
            print('Negatives: %s'%sum(negatives))

        #common to all subplots
        pl.xlabel('Time(in days starting in July)')
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
    #normal case
    testModel(vBS_oc=np.array([1.2]),subplots=[['E','A1+A2']])

    #W->0 test
    testModel(vBS_oc=np.array([1.2]),precipitations=[0 if d>500 else 15. for d in range(0,(op.end_date - op.start_date).days)],subplots=[['E','L'],['W']])

    #T->0
    time_range=range(0,(op.end_date - op.start_date).days)
    T=op.T#save the original function to be able restore it later.
    op.T=interpolate.InterpolatedUnivariateSpline(time_range,[30.*(1. - t/float(max(time_range)) )+273.15 for t in time_range ])
    testModel(vBS_oc=np.array([1.2]),subplots=[['A1+A2'],['T']])
    op.T=T#restore the original function back

    #cimsim vs spaa
    #diameter:9.5, height:5.8, type: circular, sun exposure:0.9
    #TODO:we should change the start/end dates to match the ones used in cimsim, but in this case is ok, because the container is empty anyways.
    testModel(vBS_oc=np.array([0.41]),vBS_os=np.array([math.pi* 4.75**2]),subplots=[['spaavscimsim']])

    #against Indices
    vBS_os=math.pi*np.array([42.,52.,62.])
    testModel(BS_o=0.1,vBS_oc=np.array([0.1,0.6,8.3]),vBS_os=vBS_os,vBS_od=np.array([0.0,0.0,1.0]),subplots=[['nLvsI']])
    testModel(BS_o=0.0,vBS_oc=np.array([0.1,0.6,8.3]),vBS_os=vBS_os,vBS_od=np.array([0.0,0.0,1.0]),subplots=[['nLvsI']])
    #against Ovitraps
    testModel(BS_o=0.4,vBS_oc=np.array([0.1,0.6,8.3]),vBS_os=vBS_os,vBS_od=np.array([0.1,0.5,0.4]),subplots=[['nEvsO']])

    #performace
    testModel(BS_o=0.1,vBS_oc=np.array([0.1,0.6,8.3]),vBS_os=vBS_os,vBS_od=np.array([0.0,0.0,1.0]),subplots=[['E'],['b'],['c'],['n']])

    #Different container surfaces
    '''
    for i in range(1,4):
        vBS_oc=np.array([10.0])
        op.BS_a=50.0/np.dot(op.vBS_od,op.vBS_oc)
        testModel(BS_o=1.0,vBS_oc=vBS_oc,vBS_os=math.pi*np.array([float(i)**2]),vBS_od=np.array([1.]),subplots=[['E'],['lwE'],['W'],['clc']])
    '''

    #with a various of types of containers
    '''
                      #Plant plate                  #2-l bottle in half               #dog plate                  #water tank               #'pelopincho'                  #Piscine
    vBS_os=np.array([getSurface(r=24.5/2.)        ,getSurface(r=10.15/2.)           ,getSurface(r=14./2.)       ,getSurface(r=55.)        ,getSurface(x=155.,y=107.)     ,getSurface(x=700.,y=345.)         ])
    vBS_oc=np.array([getCapacity(r=24.5/2.,z=3.084),getCapacity(r=10.15/2.,z=33.6/2.),getCapacity(r=14./2.,z=5.),getCapacity(r=55.,z=145.),getCapacity(x=155.,y=107.,z=30.),getCapacity(x=700.,y=345.,z=135.) ])
    op.ws_s=0.2
    testModel(BS_o=1.0,vBS_oc=vBS_oc,vBS_os=vBS_os,vBS_od=np.array([0.5,0.5,0.0,0.0,0.0,0.0]),subplots=[['E','A1+A2','nT','np'],['nW']])
    '''
    pl.show()
