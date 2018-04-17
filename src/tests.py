import sys
import math
import utils
import fourier
import datetime
import pylab as pl
import numpy as np
import matplotlib
import matplotlib.dates
from scipy import interpolate
from scipy.stats import stats
from config import Configuration
import otero_precipitation as op
from utils import getSurface,getCapacity#not sure if this is a good practice

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

def testRESvsOldRES(new_RES,old_RES_filename):
    old_RES_start_date=datetime.datetime.strptime(open(old_RES_filename,'r').readline().split(',')[0],'%Y-%m-%d').date()
    old_RES_end_date=datetime.datetime.strptime(open(old_RES_filename,'r').readlines()[-1].split(',')[0],'%Y-%m-%d').date()
    assert old_RES_start_date==op.start_date,'Old Result and new result must start at the same date'
    print("||old_RES-new_RES|| = "),
    old_RES=utils.loadResults(old_RES_filename,op.start_date)
    new_RES=utils.getDailyResults(op.getTimeRange(),new_RES,op.start_date,old_RES_end_date+datetime.timedelta(days=1))#the last date save is one day less than the end_date#TODO:not sure if this is correct or just a hack
    print(np.linalg.norm(old_RES-new_RES))


def testSaveLoadResults(time_range,RES):
    filename=utils.saveResults(time_range,RES,op.start_date,op.end_date)
    same_RES=utils.loadResults(filename,op.start_date)
    RES=utils.getDailyResults(time_range,RES,op.start_date,op.end_date)
    print('%s vs %s'%(len(RES),len(same_RES)))
    print(np.linalg.norm(RES-same_RES))

calls=np.array([0]*len(op.getTimeRange()))
negatives=np.array([0]*len(op.getTimeRange()))
def decoratedEquations(Y,t):
    dY=op.diff_eqs(Y,t)
    time_range=op.getTimeRange()
    calls[(np.abs(time_range-t)).argmin()]+=1
    if(np.any(Y<0)): negatives[(np.abs(time_range-t)).argmin()]+=1
    return dY

def normalize(values):#TODO:change name.
    return (values-values.min())/(values.max()-values.min())
#convenience method
def normalizeIfAsked(values,subplot):
    values=np.array(values)
    if('normalized' in subplot):
        return normalize(values)
    else:
        return values
def subData(time_range,Y,date_range,an_start_date):
    #get the index of an_start_date
    index=None
    for i,date in enumerate(date_range):
        if(date.date()==an_start_date):
            index=i
            break#conserve the first one.
    return time_range[index:],Y[index:,:],date_range[index:]

def configure(configuration):
    op.BS_a=configuration.getFloat('breeding_site','amount')
    op.vBS_oc=configuration.getArray('breeding_site','outside_capacity')#in litres
    op.vBS_ic=configuration.getArray('breeding_site','inside_capacity')#in litres
    op.vBS_od=configuration.getArray('breeding_site','outside_distribution')#distribution of BS outside # the sum of ouside and in-
    op.vBS_id=configuration.getArray('breeding_site','inside_distribution')#distribution of BS inside   #  side must be equal to 1
    op.vBS_os=configuration.getArray('breeding_site','outside_surface')#in cm^2
    op.n,op.m=len(op.vBS_od),len(op.vBS_id)
    op.ws_s=configuration.getFloat('weather','wind_shield')#wind shield in [0,1]
    #Cordoba
    op.location={'name':configuration.getString('location','name'),'station':configuration.getString('weather','station'),'zones':list(configuration.getString('location','zones'))}
    op.start_date=configuration.getDate('simulation','start_date')
    op.end_date=configuration.getDate('simulation','end_date')

    op.AEDIC_INDICES_FILENAME='data/private/Indices aedicos Historicos '+op.location['name']+'.xlsx'
    op.WEATHER_STATION_DATA_FILENAME='data/public/wunderground_'+op.location['station']+'.csv'

    precipitations = utils.getPrecipitationsFromCsv(op.WEATHER_STATION_DATA_FILENAME,op.start_date,op.end_date)
    op.p=op.getAsLambdaFunction(op.aps,precipitations)

    wind_speed=utils.getMeanWindSpeedFromCsv(op.WEATHER_STATION_DATA_FILENAME,op.start_date,op.end_date)
    op.ws=fourier.fourier(wind_speed,50)

    op.T=interpolate.InterpolatedUnivariateSpline(range(0,(op.end_date - op.start_date).days),utils.getAverageTemperaturesFromCsv(op.WEATHER_STATION_DATA_FILENAME,op.start_date,op.end_date))


def testModel(configuration, p=None,T=None,subplots=[['E','L'],['W']],plot_start_date=None):
    configure(configuration)
    if(p):
        op.p=p
    if(T):
        op.T=T
    initial_condition=configuration.getArray('simulation','initial_condition')
    time_range,INPUT,RES=op.solveEquations(initial_condition=initial_condition,equations=decoratedEquations)
    if(len(sys.argv)>1 and sys.argv[1]=='save' and p==None and T==None):#if asked save, but not with tampered p or T functions
        results_filename=utils.saveResults(time_range,RES,op.start_date,op.end_date)
        configuration.save(results_filename.replace('.csv','.cfg'))

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
        if ('E' in subplot): pl.plot(date_range,normalizeIfAsked(RES[:,op.EGG],subplot), '-k', label='E')
        if ('L' in subplot): pl.plot(date_range,normalizeIfAsked(RES[:,op.LARVAE],subplot), '-r', label='L')
        if ('P' in subplot): pl.plot(date_range,normalizeIfAsked(RES[:,op.PUPAE],subplot), '-g', label='P')
        if ('A1' in subplot): pl.plot(date_range,normalizeIfAsked(RES[:,op.ADULT1],subplot), '-b', label='A1')
        if ('A2' in subplot): pl.plot(date_range,normalizeIfAsked(RES[:,op.ADULT2],subplot), '-m', label='A2')
        if ('A1+A2' in subplot): pl.plot(date_range,normalizeIfAsked(RES[:,op.ADULT2]+RES[:,op.ADULT1],subplot), '-m', label='A1+A2')
        if ('LI' in subplot):
            ri_days, ris=utils.getIndexesForPlot(op.AEDIC_INDICES_FILENAME,op.start_date,0,5)
            pl.plot([datetime.timedelta(days=d)+datetime.datetime.combine(op.start_date,datetime.time()) for d in ri_days], normalizeIfAsked(ris,subplot), '^y', label='Recip. Indices',clip_on=False, zorder=100,markersize=8)

        if('O' in subplot):
            for i in subplot['O']:
                ovitrap_eggs=utils.getOvitrapEggsFromCsv('data/private/Datos sensores de oviposicion.NO.csv',op.start_date,op.end_date,i)
                if('normalized' in subplot):ovitrap_eggs=np.array([e/max(ovitrap_eggs) if e!=None else None for e in ovitrap_eggs])#since we have None values normalized won't work
                pl.plot([datetime.timedelta(days=d)+datetime.datetime.combine(op.start_date,datetime.time()) for d in range(0,len(ovitrap_eggs))], ovitrap_eggs, '^', label='Ovitrap %s eggs'%i,clip_on=False, zorder=100,markersize=8)

        if('lwE' in subplot):
            lwE=np.array([RES[(np.abs(time_range-t)).argmin(),op.EGG]-RES[(np.abs(time_range-(t-7))).argmin(),op.EGG] for t in time_range])
            if('normalized' in subplot):#not same normalize as a
                    lwE[lwE<0]=0.#replace negatives with zeros
            pl.plot(date_range, normalizeIfAsked(lwE,subplot), '-', label='E(t)-E(t-7)')
        pl.ylabel('')

        #Complete lifecycle
        if('clc' in subplot):
            pl.plot(date_range,[sum([1./op.R_D(stage,op.T(t)) for stage in [op.EGG,op.LARVAE,op.PUPAE,op.ADULT1,op.ADULT2]]) for  t in time_range],label='Complete life cicle(from being an egg to the second oviposition)')
            pl.ylabel('days')
        #Water in containers(in L)
        if ('W' in subplot):
            for i in range(0,op.n):
                pl.plot(date_range,normalizeIfAsked(RES[:,op.WATER+i],subplot), label='W(t) for %sL, %scm^2, %s%%'%(op.vBS_oc[i],op.vBS_os[i],op.vBS_od[i]*100.) )
            pl.ylabel('Litres')

        #spaa vs cimsim
        if ('spaavscimsim' in subplot):
            for i in range(0,op.n):
                pl.plot(time_range,RES[:,op.WATER+i]*1000.0/op.vBS_os[i], label='W(t) for %sL, %scm^2, %s%%'%(op.vBS_oc[i],op.vBS_os[i],op.vBS_od[i]*100.) )#L->ml->mm->cm
            pl.plot(utils.getValuesFromCsv('data/test/cimsim_containers_2015_se09.csv',op.start_date,op.end_date,1,verbose=False),label='CIMSiM')

        #Temperature in K
        if ('T' in subplot):
            pl.plot(date_range,normalizeIfAsked([op.T(t) for t in time_range],subplot), label='Temperature')
            pl.ylabel('K')

        #precipitations(in mm.)
        if ('p' in subplot):
            pl.plot(date_range,normalizeIfAsked([op.p(t) for t in time_range],subplot),'-b', label='p(t)')
            pl.ylabel('mm./day')

        #Wind Speed(in km/h.)
        if ('ws' in subplot):
            pl.plot(date_range,normalizeIfAsked([op.ws(t) for t in time_range],subplot), label='ws(t)')
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
    print('Larvae Survival to desiccation:')
    dry_time=times_dry_container[0]
    print(' W<epsilon at t=%f with epsilon=%f'%(dry_time,epsilon) )
    print(' L<epsilon at t=%f with epsilon=%f'%(dry_time,epsilon) )
    print(' L=%s at t=%f'%(L[np.where(time_range==dry_time)],dry_time))

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
    if(len(sys.argv)>1 and sys.argv[1]=='compare'):
        filenames=[]
        if(len(sys.argv)>2):
            filenames=sys.argv[2:]#command line
        else:
            filenames=sys.stdin.read().split('\n')[:-1]#command pipe

        if(len(filenames)>0):
            for filename in filenames:
                configuration=Configuration(filename.replace('.csv','.cfg'))
                configure(configuration)
                initial_condition=configuration.getArray('simulation','initial_condition')
                time_range,INPUT,RES=op.solveEquations(initial_condition=initial_condition)
                testRESvsOldRES(RES,filename)
            quit()

    config=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'outside_capacity':[1.2],
            'outside_surface':math.pi*np.array([5.25**2]),
            'outside_distribution':[1.],
            'inside_distribution':[0]
            },
        'simulation':{
            'initial_condition':[100.0, 0.0,0.0,0.0,0.0]+ [0]
            }
        })
    #normal case
    testModel(config,subplots=[['E','A1+A2','normalized']])

    #W->0 test
    precipitations=[0 if d>500 else 15. for d in range(0,(op.end_date - op.start_date).days)]
    p=op.getAsLambdaFunction(op.aps,precipitations)
    testModel(config,p=p,subplots=[['E','L','normalized'],['W']])

    #T->0
    time_range=range(0,(config.getDate('simulation','end_date') - config.getDate('simulation','start_date')).days )
    T=interpolate.InterpolatedUnivariateSpline(time_range,[30.*(1. - t/float(max(time_range)) )+273.15 for t in time_range ])
    testModel(config,subplots=[['A1+A2'],['T']],T=T)

    #cimsim vs spaa
    config=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'outside_capacity':[0.41],
            'outside_surface':math.pi*np.array([4.75**2]),
            'outside_distribution':[1.0],
            'inside_distribution':[0]
            },
        'simulation':{
            'initial_condition':[100.0, 0.0,0.0,0.0,0.0]+ [0]
            }
        })
    #diameter:9.5, height:5.8, type: circular, sun exposure:0.9
    #TODO:we should change the start/end dates to match the ones used in cimsim, but in this case is ok, because the container is empty anyways.
    testModel(config,subplots=[['spaavscimsim']])

    #against Indices
    vBS_os=math.pi*np.array([42.,52.,62.])
    config=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'outside_capacity':[0.1,0.6,8.3],
            'outside_surface':vBS_os,
            'outside_distribution':[0,0,0.1],
            'inside_distribution':[0.9]
            },
        'simulation':{
            'initial_condition':[100.0, 0.0,0.0,0.0,0.0]+ [0 for x in vBS_os]
            }
        })
    testModel(config,subplots=[['L','LI','normalized']])

    config=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'outside_capacity':[0.1,0.6,8.3],
            'outside_surface':vBS_os,
            'outside_distribution':[0,0,0.0],
            'inside_distribution':[1.0]
            },
        'simulation':{
            'initial_condition':[100.0, 0.0,0.0,0.0,0.0]+ [0 for x in vBS_os]
            }
        })
    testModel(config,subplots=[['L','LI','normalized']])

    #performace
    testModel(config,subplots=[['E'],['b'],['c'],['n']])


    #with a various of types of containers

                      #Plant plate                  #2-l bottle in half               #dog plate                  #water tank               #'pelopincho'                  #Piscine
    vBS_os=np.array([getSurface(r=24.5/2.)        ,getSurface(r=10.15/2.)           ,getSurface(r=14./2.)       ,getSurface(r=55.)        ,getSurface(x=155.,y=107.)     ,getSurface(x=700.,y=345.)         ])
    vBS_oc=np.array([getCapacity(r=24.5/2.,z=3.084),getCapacity(r=10.15/2.,z=33.6/2.),getCapacity(r=14./2.,z=5.),getCapacity(r=55.,z=145.),getCapacity(x=155.,y=107.,z=30.),getCapacity(x=700.,y=345.,z=135.) ])

    config=Configuration('resources/otero_precipitation.cfg',{
        'breeding_site':{
            'outside_capacity':vBS_oc,
            'outside_surface':vBS_os,
            'outside_distribution':[0.5,0.5,0.0,0.0,0.0,0.0],
            'inside_distribution':[0]
            },
        'weather':{
            'wind_shield':0.2
        },
        'simulation':{
            'initial_condition':[100.0, 0.0,0.0,0.0,0.0]+ [0 for x in vBS_os]
        }
    })
    testModel(config,subplots=[['E','A1+A2','T','p','normalized'],['W']],plot_start_date=datetime.date(2018,1,1))

    #*****9*****
    #ovitrap:9 pid:2382 od:[ 0.03088072  0.20904943  0.23383199  0.16713309  0.17310652  0.11768087] id:[ 0.06831738] ws_s:0.031265688907 Error:0.0765284863715 len:11.0 Error/len: 0.00695713512468
    vBS_od=np.array([0.03088072,0.20904943,0.23383199,0.16713309,0.17310652,0.11768087])
    BS_o=np.sum(vBS_od)
    print(BS_o)
    config=Configuration('resources/otero_precipitation.cfg',{
        'breeding_site':{
            'outside_capacity':vBS_oc,
            'outside_surface':vBS_os,
            'outside_distribution':vBS_od,
            'inside_distribution':[0.06831738]
            },
        'weather':{
            'wind_shield':0.031265688907
        },
        'simulation':{
            'start_date':datetime.date(2017,7,1),
            'end_date':datetime.date(2018,4,5),
            'initial_condition':[100.0, 0.0,0.0,0.0,0.0]+ [0 for x in vBS_os]
        }
    })
    testModel(config,subplots=[['E','P','A1+A2','normalized'],{'lwE':'','O':[9],'normalized':''}])

    #*****4*****
    #ovitrap:4 pid:18743 od:[ 0.18299322  0.20899391  0.07332913  0.15454651  0.14291156  0.0308964 ] id:[ 0.20632926] ws_s:0.491606121558 BS_a:2594.27715109 Error:34425.9670772 len:18.0 Error/len: 1912.553
    vBS_od=np.array([0.18299322,0.20899391,0.07332913,0.15454651,0.14291156,0.0308964])
    BS_o=np.sum(vBS_od)
    print(BS_o)
    config=Configuration('resources/otero_precipitation.cfg',{
        'breeding_site':{
            'amount':2595,
            'outside_capacity':vBS_oc,
            'outside_surface':vBS_os,
            'outside_distribution':vBS_od,
            'inside_distribution':[0.20632927]
            },
        'weather':{
            'wind_shield':0.421606121558
        },
        'simulation':{
            'start_date':datetime.date(2017,7,1),
            'end_date':datetime.date(2018,4,5),
            'initial_condition':[100.0, 0.0,0.0,0.0,0.0]+ [0 for x in vBS_os]
        }
    })
    testModel(config,subplots=[['E','P','A1+A2','normalized'],{'lwE':'','O':[4],'normalized':''}])

    #*****3*****
    #ovitrap:3 pid:18743 od:[ 0.07533379  0.35492456  0.0164825   0.04007676  0.08755963  0.0680057 ] id:[ 0.35761705] ws_s:0.895738915951 BS_a:3132.19610422 Error:5057.73452148 len:20.0 Error/len: 252.886726074
    vBS_od=np.array([0.07533379 , 0.35492456 , 0.0164825   ,0.04007676 , 0.08755963 , 0.0680057])
    BS_o=np.sum(vBS_od)
    print(BS_o)
    config=Configuration('resources/otero_precipitation.cfg',{
        'breeding_site':{
            'amount':3132,
            'outside_capacity':vBS_oc,
            'outside_surface':vBS_os,
            'outside_distribution':vBS_od,
            'inside_distribution':[0.35761706]
            },
        'weather':{
            'wind_shield':0.31738915951
        },
        'simulation':{
            'start_date':datetime.date(2017,7,1),
            'end_date':datetime.date(2018,4,5),
            'initial_condition':[100.0, 0.0,0.0,0.0,0.0]+ [0 for x in vBS_os]
        }
    })
    testModel(config,subplots=[['E','P','A1+A2','normalized'],{'lwE':'','O':[3],'normalized':''}])

    #*****4 but just to compare with something*****
    vBS_od=np.array([ 0.01 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ])
    BS_o=np.sum(vBS_od)
    print(BS_o)
    config=Configuration('resources/otero_precipitation.cfg',{
        'breeding_site':{
            'amount':2595,
            'outside_capacity':vBS_oc,
            'outside_surface':vBS_os,
            'outside_distribution':vBS_od,
            'inside_distribution':[0.99]
            },
        'weather':{
            'wind_shield':0.421606121558
        },
        'simulation':{
            'start_date':datetime.date(2017,7,1),
            'end_date':datetime.date(2018,4,5),
            'initial_condition':[100.0, 0.0,0.0,0.0,0.0]+ [0 for x in vBS_os]
        }
    })
    testModel(config,subplots=[['E','P','A1+A2','normalized'],{'lwE':'','O':[4],'normalized':''}])



    pl.show()
