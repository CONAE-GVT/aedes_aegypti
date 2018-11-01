#coding: utf-8
import os
import re
import sys
import math
import utils
import datetime
import numpy as np
from scipy import interpolate
from scipy.stats import stats
from config import Configuration
from otero_precipitation import Model
from utils import getSurface,getCapacity#not sure if this is a good practice
from equations import diff_eqs

def printCorrelation():
    time_range,INPUT,RES=model.solveEquations()
    vT=np.vectorize(model.T)
    print('(Pearson s correlation coefficient,2-tailed p-value): ' + str(stats.pearsonr(RES[:,LARVAE], vT(time_range))) )

def getLarvaeSurvivalDry(time_range,vW,L,time_window):
    #print(W[:,0]<1e-3)
    #print(time_range[W<1e-3])
    epsilon=1e-4
    times_dry_container=[t for i,t in enumerate(time_range) if np.dot(model.vBS_od,vW[i,:])<epsilon and time_window['start']<t<time_window['end'] ]
    times_death_larvae= [t for i,t in enumerate(time_range) if L[i]<epsilon and time_window['start']<t<time_window['end']]
    if not times_dry_container or not times_death_larvae:
        return
    print('Larvae Survival to desiccation:')
    dry_time=times_dry_container[0]
    print(' W<epsilon at t=%f with epsilon=%f'%(dry_time,epsilon) )
    print(' L<epsilon at t=%f with epsilon=%f'%(dry_time,epsilon) )
    print(' L=%s at t=%f'%(L[np.where(time_range==dry_time)],dry_time))

def plotTimeWaterEggs(configuration):#from https://jakevdp.github.io/PythonDataScienceHandbook/04.12-three-dimensional-plotting.html
    from scipy import interpolate
    model=Model(configuration)
    time_range,INPUT,RES=model.solveEquations(equations=diff_eqs)
    WATER=model.parameters.WATER
    eggs=interpolate.SmoothBivariateSpline(time_range,RES[:,WATER[0]], RES[:,0])
    xline = time_range
    yline = RES[:,WATER[0]]
    zline = np.array([eggs(time_range[i],RES[i,WATER[0]]) for i,t in enumerate(time_range)]).flatten()
    utils.plot3D(xline,yline,zline)
    utils.showPlot()


def testSaveLoadResults(model):
    filename=model.save()
    same_RES=utils.loadResults(filename,model.start_date)
    RES=utils.getDailyResults(model.time_range,model.Y,model.start_date,model.end_date)
    print('%s vs %s'%(len(RES),len(same_RES)))
    print(np.linalg.norm(RES-same_RES))

class DecoratedEquations:
    def __init__(self,model,diff_eqs):
        self.model=model
        self.diff_eqs=diff_eqs

    def __call__(self,Y,t,parameters):
        dY=self.diff_eqs(Y,t,parameters)
        time_range=self.model.time_range
        #account calls
        if(parameters.calls is None):
            parameters.calls=np.array([0]*len(time_range))
        parameters.calls[(np.abs(time_range-t)).argmin()]+=1
        #account negatives
        if(parameters.negatives is None):
            parameters.negatives=np.array([0]*len(time_range))
        if(np.any(Y<0)): parameters.negatives[(np.abs(time_range-t)).argmin()]+=1
        return dY

def testModel(configuration, p=None,T=None,subplots=[['E','L'],['W']],plot_start_date=None):
    model=Model(configuration)

    if(p):
        model.parameters.weather.p=p
    if(T):
        model.parameters.weather.T=T

    model.parameters.calls=None
    model.parameters.negatives=None
    time_range,INPUT,RES=model.solveEquations(equations=DecoratedEquations(model,diff_eqs) )
    if(len(sys.argv)>1 and sys.argv[1]=='save' and p==None and T==None):#if asked save, but not with tampered p or T functions
        model.save()

    utils.plot(model,subplots,plot_start_date,title=configuration.getString('location','name'))


def compare(old_RES_filename,model):
    old_RES_start_date=datetime.datetime.strptime(open(old_RES_filename,'r').readline().split(',')[0],'%Y-%m-%d').date()
    old_RES_end_date=datetime.datetime.strptime(open(old_RES_filename,'r').readlines()[-1].split(',')[0],'%Y-%m-%d').date()
    assert old_RES_start_date==model.start_date,'Old Result and new result must start at the same date'
    old_RES=utils.loadResults(old_RES_filename,model.start_date)
    time_range,INPUT,new_RES=model.solveEquations(equations=diff_eqs)
    new_RES=utils.getDailyResults(time_range,new_RES,model.start_date,old_RES_end_date+datetime.timedelta(days=1))#the last date save is one day less than the end_date#TODO:not sure if this is correct or just a hack
    print('||old_RES-new_RES|| = %s'% np.linalg.norm(old_RES-new_RES) )

def runComparison():
    filenames=[]
    if(len(sys.argv)>2):
        filenames=sys.argv[2:]#command line
    else:
        filenames=sys.stdin.read().split('\n')[:-1]#command pipe

    if(len(filenames)>0):
        for filename in filenames:
            configuration=Configuration(filename.replace('.csv','.cfg'))
            model=Model(configuration)
            compare(filename,model)
        quit()

def runTestCases():

    config=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'height':[13.86],
            'surface':math.pi*np.array([5.25**2]),
            'distribution':[1.],
            'manually_filled':[0.5]
            },
        'simulation':{
            'initial_condition':[100.]*1 + [0.]*1 +[0.]*1 + [0.,0.]+ [9]
            },
        'biology':{
            'alpha0':[1.5]
            }
        })
    #normal case
    testModel(config,subplots=[['E','A1+A2',[utils.safeAdd,utils.normalize] ]])

    #W->0 test
    tmp=Model()#kind of a hack
    precipitations=[0 if d>500 else 15. for d in range(0,(tmp.end_date - tmp.start_date).days)]
    p=tmp.parameters.weather.getAsLambdaFunction(tmp.parameters.weather.aps,precipitations)
    testModel(config,p=p,subplots=[['E','L',[utils.safeAdd,utils.normalize] ],['W']])

    #T->0
    time_range=range(0,(config.getDate('simulation','end_date') - config.getDate('simulation','start_date')).days )
    T=interpolate.InterpolatedUnivariateSpline(time_range,[30.*(1. - t/float(max(time_range)) )+273.15 for t in time_range ])
    testModel(config,subplots=[['A1+A2'],['T']],T=T)

    #cimsim vs spaa
    config=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'height':[5.8],
            'surface':math.pi*np.array([4.75**2]),
            'distribution':[1.0],
            'manually_filled':[0.5]
            },
        'simulation':{
            'start_date':datetime.date(2015,7,15),
            'end_date': datetime.date(2015,12,31),
            'initial_condition':[100.]*1 + [0.]*1 +[0.]*1 + [0.,0.]+ [0]
            },
        'biology':{
            'alpha0':[1.5]
            }
        })
    #diameter:9.5, height:5.8, type: circular, sun exposure:0.9
    #TODO:we should change the start/end dates to match the ones used in cimsim, but in this case is ok, because the container is empty anyways.
    testModel(config,subplots=[['spaavscimsim']])

    #against Indices
    vBS_s=math.pi*np.array([50, 42.,52.,62.])
    config=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'height':[10,9.1],
            'surface':vBS_s[0:2],
            'distribution':[0.9,0.1],
            'manually_filled':[0.5,0]
            },
        'simulation':{
            'initial_condition':[100.]*2 + [0.]*2 +[0.]*2 + [0.,0.]+ [0 for x in vBS_s[0:2]]
            },
        'biology':{
            'alpha0':[1.5]*2
            }
        })
    testModel(config,subplots=[['L','LI',[utils.safeAdd,utils.normalize] ]])

    config=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'height':[10],
            'surface':[50],
            'distribution':[1.0],
            'manually_filled':[0.5]
            },
        'simulation':{
            'initial_condition':[100.]*1 + [0.]*1 +[0.]*1 + [0.,0.] + [9.]
            },
        'biology':{
            'alpha0':[1.5]*1
            }
        })
    testModel(config,subplots=[['L','LI',[utils.safeAdd,utils.normalize] ]])

    #performace
    testModel(config,subplots=[['E'],['b'],['c'],['n']])


    #with a various of types of containers

                      #Plant plate                  #2-l bottle in half               #dog plate                  #water tank               #'pelopincho'                  #Piscine
    vBS_s=np.array([50, getSurface(r=24.5/2.)        ,getSurface(r=10.15/2.)           ,getSurface(r=14./2.)       ,getSurface(r=55.)        ,getSurface(x=155.,y=107.)     ,getSurface(x=700.,y=345.)         ])
    vBS_h=np.array([10, 3.084                        ,33.6/2.                          ,5.                         ,145                      ,30                            ,135.])

    config=Configuration('resources/otero_precipitation.cfg',{
        'breeding_site':{
            'height':vBS_h[1:3],
            'surface':vBS_s[1:3],
            'distribution':[0.5,0.5],
            'manually_filled':[0.5,0.]
            },
        'simulation':{
            'initial_condition':[100.]*2 + [0.]*2 +[0.]*2 + [0.,0.]+ [0]*len(vBS_s[0:2])
        },
        'biology':{
            'alpha0':[1.5]*2
        }
    })
    testModel(config,subplots=[['E','A1+A2','T','p','RH',[utils.safeAdd,utils.normalize] ],['W'],['mf']],plot_start_date=datetime.date(2018,1,1))

    #*****9*****
    #ovitrap:9 pid:2382 od:[ 0.03088072  0.20904943  0.23383199  0.16713309  0.17310652  0.11768087] id:[ 0.06831738] ws_s:0.031265688907 Error:0.0765284863715 len:11.0 Error/len: 0.00695713512468
    vBS_d=np.array([0.06831738, 0.03088072,0.20904943,0.23383199,0.16713309,0.17310652,0.11768087])
    config=Configuration('resources/otero_precipitation.cfg',{
        'breeding_site':{
            'height':vBS_h,
            'surface':vBS_s,
            'distribution':vBS_d
            },
        'simulation':{
            'start_date':datetime.date(2017,7,1),
            'end_date':datetime.date(2018,4,5),
            'initial_condition':[100.]*7 + [0.]*7 +[0.]*7 + [0.,0.]+ [0]*len(vBS_s)
        }
    })
    testModel(config,subplots=[['E','P','A1+A2',[utils.safeAdd,utils.normalize] ],['dE','dL','dP','dA1','dA2',[utils.safeAdd]],{'lwE':'','O':[9],'f' :[utils.replaceNegativesWithZeros,utils.safeAdd,utils.safeNormalize]}])

    #*****4*****
    #ovitrap:4 pid:18743 od:[ 0.18299322  0.20899391  0.07332913  0.15454651  0.14291156  0.0308964 ] id:[ 0.20632926] ws_s:0.491606121558 BS_a:2594.27715109 Error:34425.9670772 len:18.0 Error/len: 1912.553
    vBS_d=np.array([0.20632927, 0.18299322,0.20899391,0.07332913,0.15454651,0.14291156,0.0308964])
    config=Configuration('resources/otero_precipitation.cfg',{
        'breeding_site':{
            'amount':2595,
            'height':vBS_h,
            'surface':vBS_s,
            'distribution':vBS_d
            },
        'simulation':{
            'start_date':datetime.date(2017,7,1),
            'end_date':datetime.date(2018,4,5),
            'initial_condition':[100.]*7 + [0.]*7 +[0.]*7 + [0.,0.]+ [0]*len(vBS_s)
        }
    })
    testModel(config,subplots=[['E','P','A1+A2',[utils.safeAdd,utils.normalize] ],{'lwE':'','O':[4],'f':[utils.replaceNegativesWithZeros,utils.safeAdd,utils.safeNormalize]}])

    #*****3*****
    #ovitrap:3 pid:18743 od:[ 0.07533379  0.35492456  0.0164825   0.04007676  0.08755963  0.0680057 ] id:[ 0.35761705] ws_s:0.895738915951 BS_a:3132.19610422 Error:5057.73452148 len:20.0 Error/len: 252.886726074
    vBS_d=np.array([0.35761706, 0.07533379 , 0.35492456 , 0.0164825   ,0.04007676 , 0.08755963 , 0.0680057])
    config=Configuration('resources/otero_precipitation.cfg',{
        'breeding_site':{
            'amount':3132,
            'height':vBS_h,
            'surface':vBS_s,
            'distribution':vBS_d
            },
        'simulation':{
            'start_date':datetime.date(2017,7,1),
            'end_date':datetime.date(2018,4,5),
            'initial_condition':[100.]*7 + [0.]*7 +[0.]*7 + [0.,0.]+ [0]*len(vBS_s)
        }
    })
    testModel(config,subplots=[['E','P','A1+A2',[utils.safeAdd,utils.normalize] ],{'lwE':'','O':[3],'f':[utils.replaceNegativesWithZeros,utils.safeAdd,utils.safeNormalize]}])

    #*****4 but just to compare with something*****
    vBS_d=np.array([0.99, 0.01])
    config=Configuration('resources/otero_precipitation.cfg',{
        'breeding_site':{
            'amount':2595,
            'height':vBS_h[0:2],
            'surface':vBS_s[0:2],
            'distribution':vBS_d[0:2],
            'manually_filled':[0.5,0.]
            },
        'simulation':{
            'start_date':datetime.date(2017,7,1),
            'end_date':datetime.date(2018,4,5),
            'initial_condition':[100.]*2 + [0.]*2 +[0.]*2 + [0.,0.]+ [0]*len(vBS_s[0:2])
        },
        'biology':{
            'alpha0':[1.5]*2
        }
    })
    testModel(config,subplots=[['E','P','A1+A2',[utils.safeAdd,utils.normalize] ],{'lwE':'','O':[4],'f':[utils.replaceNegativesWithZeros,utils.safeAdd,utils.safeNormalize]}])

    #plotTimeWaterEggs(config)
    utils.showPlot()

def getEquilibriumL0(u,T_t):
    from equations import vR_D
    elr,lpr,par,ovr1,ovr2=vR_D(T_t)
    egn=63.0
    me=0.01#mortality of the egg, for T in [278,303]
    ml=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#mortality of the larvae, for T in [278,303]
    mp=0.01 + 0.9725 * math.exp(-(T_t-278.0)/2.7035)#death of pupae
    ef=0.83#emergence factor
    ma=0.09#for T in [278,303]
    alpha0=1.5
    BS=1.
    alpha=alpha0/BS
    return 1./alpha *( (elr * u * egn*(ma*ovr1+ovr2*ovr1)*par*ef*lpr)/(2*ma*(me+elr*u)*(ma+ovr1)*(mp+par) ) - (ml+lpr) )

def testMacia():#MACIÁ (2009)
    configuration=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'amount':1,
            'outside_capacity':[],
            'outside_surface':[],
            'outside_distribution':[],
            'inside_distribution':[1.],
            'inside_capacity':[0.175]
            },
        'simulation':{
            'start_date':datetime.date(2017,3,1),
            'end_date':datetime.date(2017,4,15),
            'initial_condition':[0.]*1 + [32.]*1 +[0.]*1 + [0.,0.]
            },
        'biology':{
            'alpha0':[1.5]
            }
        })

    model=Model(configuration)
    #lab environment
    T=lambda t: 26.+273.15#26 +/- 2
    model.parameters.weather.T=T
    time_range,INPUT,RES=model.solveEquations()
    LARVAE,PUPAE,WATER=model.parameters.LARVAE,model.parameters.PUPAE,model.parameters.WATER
    vBS_ic=model.parameters.vBS_ic
    L_initial,P_final,D_initial=RES[0,LARVAE],RES[-1,PUPAE],RES[0,LARVAE]/np.concatenate((RES[0,WATER],vBS_ic) )
    print('D_initial: %s, L_initial:%s, P_final:%s, survival:%s'%(D_initial,L_initial,P_final,P_final/L_initial))

def runOviShow(folder):
    config_filenames=[filename for filename in os.listdir(folder) if filename.endswith('.cfg')]
    for config_filename in config_filenames:
        for ovitrap_id in re.findall(r'.*_ovi([0-9]+)\.cfg',config_filename):
            testModel(Configuration(folder+'/'+config_filename),subplots=[ ['E','P','A1+A2',[utils.safeAdd,utils.normalize] ],{'lwE':'','O':[int(ovitrap_id)+1],'f':[utils.replaceNegativesWithZeros,utils.safeAdd,utils.safeNormalize]}])

    utils.showPlot()

def runLab():
    testMacia()

import pylab as pl
def compareWeather(gfs_weather_filename,cfs_weather_filename):
        start_date=datetime.date.today()
        end_date=datetime.date.today()+datetime.timedelta(7)
        temperatures_gfs= np.array(utils.getAverageTemperaturesFromCsv(gfs_weather_filename,start_date,end_date))
        temperatures_cfs= np.array(utils.getAverageTemperaturesFromCsv(cfs_weather_filename,start_date,end_date))
        print(gfs_weather_filename.split('/')[-1].replace('.csv','') + ':')

        delta_temperature=temperatures_gfs-temperatures_cfs
        print('T: %s +/- %s'%(np.mean(delta_temperature),np.std(delta_temperature)) )

        #pl.figure()
        RH_gfs= np.array(utils.getRelativeHumidityFromCsv(gfs_weather_filename,start_date,end_date))
        RH_cfs= np.array(utils.getRelativeHumidityFromCsv(cfs_weather_filename,start_date,end_date))

        delta_RH=RH_gfs-RH_cfs
        print('RH:  %s +/- %s'%(np.mean(delta_RH),np.std(delta_RH)) )


def runProject():
    config=Configuration('resources/otero_precipitation.cfg',
        {'simulation':{
            'start_date':datetime.date(2017,7,1),
            'end_date':datetime.date.today()+datetime.timedelta(6),
        }
        })
    for location in utils.getLocations():
        config.config_parser.set('location','name',location+'.full')
        testModel(config,subplots=[['E','A1+A2',[utils.safeAdd,utils.normalize] ]])

    utils.showPlot()

def createMaps():
    start_date=datetime.date(2017,7,1)
    end_date=datetime.date.today()+datetime.timedelta(6)
    config=Configuration('resources/otero_precipitation.cfg',
        {'simulation':{
            'start_date':start_date,
            'end_date':end_date,
        }
        })
    for location in utils.getLocations():
        config.config_parser.set('location','name',location+'.full')
        model=Model(config)
        dates=[datetime.date(2017,12,1),datetime.date(2018,3,1),datetime.date.today()]
        utils.createMap(model,dates)

def generateCSV(start_date,end_date):
    filename='data/public/cordoba.csv'
    temperatures=np.array(utils.getAverageTemperaturesFromCsv(filename,start_date,end_date))
    relative_humiditys=np.array(utils.getRelativeHumidityFromCsv(filename,start_date,end_date))
    precipitations = np.array(utils.getPrecipitationsFromCsv(filename,start_date,end_date))

    temperatures=np.mean(temperatures.reshape(-1,7),axis=1)
    relative_humiditys=np.mean(relative_humiditys.reshape(-1,7),axis=1)
    precipitations=np.sum(precipitations.reshape(-1,7),axis=1)#should be accumulated?

    print('Temperatura, Humedad Relativa, Precipitacion')
    for i in range(len(precipitations)):
        print('%s, %s, %s'%(temperatures[i],relative_humiditys[i],precipitations[i]))


if(__name__ == '__main__'):
    if(len(sys.argv)>1 and sys.argv[1]=='compare'):
        runComparison()
    elif(len(sys.argv)>1 and sys.argv[1]=='lab'):
        runLab()
    elif(len(sys.argv)>2 and sys.argv[1]=='show'):
        runOviShow(sys.argv[2])
    elif(len(sys.argv)>1 and sys.argv[1]=='weather'):
        gfs_folder,cfs_folder=sys.argv[2:4]
        for location in utils.getLocations():
            compareWeather(gfs_folder+'/'+location+'.forecast.csv',cfs_folder+'/'+location+'.forecast.csv')
    elif(len(sys.argv)>1 and sys.argv[1]=='project'):
        runProject()
    elif(len(sys.argv)>1 and sys.argv[1]=='map'):
        createMaps()
    elif(len(sys.argv)>1 and sys.argv[1]=='csv'):
        FORMAT='%Y-%m-%d'
        start_date,end_date= datetime.datetime.strptime(sys.argv[2],FORMAT).date(),datetime.datetime.strptime(sys.argv[3],FORMAT).date()
        assert (end_date-start_date).days%7==0,(end_date-start_date).days%7
        generateCSV(start_date,end_date)
    else:
        runTestCases()
