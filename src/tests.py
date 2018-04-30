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

def compare(new_RES,old_RES_filename,model):
    old_RES_start_date=datetime.datetime.strptime(open(old_RES_filename,'r').readline().split(',')[0],'%Y-%m-%d').date()
    old_RES_end_date=datetime.datetime.strptime(open(old_RES_filename,'r').readlines()[-1].split(',')[0],'%Y-%m-%d').date()
    assert old_RES_start_date==model.start_date,'Old Result and new result must start at the same date'
    old_RES=utils.loadResults(old_RES_filename,model.start_date)
    new_RES=utils.getDailyResults(model.time_range,new_RES,model.start_date,old_RES_end_date+datetime.timedelta(days=1))#the last date save is one day less than the end_date#TODO:not sure if this is correct or just a hack
    print('||old_RES-new_RES|| = %s'% np.linalg.norm(old_RES-new_RES) )


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

    utils.plot(model,subplots,plot_start_date)

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
            time_range,INPUT,RES=model.solveEquations()
            compare(RES,filename,model)
        quit()

def runTestCases():

    config=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'outside_capacity':[1.2],
            'outside_surface':math.pi*np.array([5.25**2]),
            'outside_distribution':[1.],
            'inside_distribution':[0]
            },
        'simulation':{
            'initial_condition':[100.]*2 + [0.]*2 +[0.]*2 + [0.,0.]+ [0]
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
            'outside_capacity':[0.41],
            'outside_surface':math.pi*np.array([4.75**2]),
            'outside_distribution':[1.0],
            'inside_distribution':[0]
            },
        'simulation':{
            'initial_condition':[100.]*2 + [0.]*2 +[0.]*2 + [0.,0.]+ [0]
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
            'initial_condition':[100.]*4 + [0.]*4 +[0.]*4 + [0.,0.]+ [0 for x in vBS_os]
            }
        })
    testModel(config,subplots=[['L','LI',[utils.safeAdd,utils.normalize] ]])

    config=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'outside_capacity':[0.1,0.6,8.3],
            'outside_surface':vBS_os,
            'outside_distribution':[0,0,0.0],
            'inside_distribution':[1.0]
            },
        'simulation':{
            'initial_condition':[100.]*4 + [0.]*4 +[0.]*4 + [0.,0.]+ [0 for x in vBS_os]
            }
        })
    testModel(config,subplots=[['L','LI',[utils.safeAdd,utils.normalize] ]])

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
            'initial_condition':[100.]*7 + [0.]*7 +[0.]*7 + [0.,0.]+ [0 for x in vBS_os]
        }
    })
    testModel(config,subplots=[['E','A1+A2','T','p',[utils.safeAdd,utils.normalize] ],['W']],plot_start_date=datetime.date(2018,1,1))

    #*****9*****
    #ovitrap:9 pid:2382 od:[ 0.03088072  0.20904943  0.23383199  0.16713309  0.17310652  0.11768087] id:[ 0.06831738] ws_s:0.031265688907 Error:0.0765284863715 len:11.0 Error/len: 0.00695713512468
    vBS_od=np.array([0.03088072,0.20904943,0.23383199,0.16713309,0.17310652,0.11768087])
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
            'initial_condition':[100.]*7 + [0.]*7 +[0.]*7 + [0.,0.]+ [0 for x in vBS_os]
        }
    })
    testModel(config,subplots=[['E','P','A1+A2',[utils.safeAdd,utils.normalize] ],['dY'],{'lwE':'','O':[9],'f' :[utils.replaceNegativesWithZeros,utils.safeAdd,utils.safeNormalize]}])

    #*****4*****
    #ovitrap:4 pid:18743 od:[ 0.18299322  0.20899391  0.07332913  0.15454651  0.14291156  0.0308964 ] id:[ 0.20632926] ws_s:0.491606121558 BS_a:2594.27715109 Error:34425.9670772 len:18.0 Error/len: 1912.553
    vBS_od=np.array([0.18299322,0.20899391,0.07332913,0.15454651,0.14291156,0.0308964])
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
            'initial_condition':[100.]*7 + [0.]*7 +[0.]*7 + [0.,0.]+ [0 for x in vBS_os]
        }
    })
    testModel(config,subplots=[['E','P','A1+A2',[utils.safeAdd,utils.normalize] ],{'lwE':'','O':[4],'f':[utils.replaceNegativesWithZeros,utils.safeAdd,utils.safeNormalize]}])

    #*****3*****
    #ovitrap:3 pid:18743 od:[ 0.07533379  0.35492456  0.0164825   0.04007676  0.08755963  0.0680057 ] id:[ 0.35761705] ws_s:0.895738915951 BS_a:3132.19610422 Error:5057.73452148 len:20.0 Error/len: 252.886726074
    vBS_od=np.array([0.07533379 , 0.35492456 , 0.0164825   ,0.04007676 , 0.08755963 , 0.0680057])
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
            'initial_condition':[100.]*7 + [0.]*7 +[0.]*7 + [0.,0.]+ [0 for x in vBS_os]
        }
    })
    testModel(config,subplots=[['E','P','A1+A2',[utils.safeAdd,utils.normalize] ],{'lwE':'','O':[3],'f':[utils.replaceNegativesWithZeros,utils.safeAdd,utils.safeNormalize]}])

    #*****4 but just to compare with something*****
    vBS_od=np.array([ 0.01 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ])
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
            'initial_condition':[100.]*7 + [0.]*7 +[0.]*7 + [0.,0.]+ [0 for x in vBS_os]
        }
    })
    testModel(config,subplots=[['E','P','A1+A2',[utils.safeAdd,utils.normalize] ],{'lwE':'','O':[4],'f':[utils.replaceNegativesWithZeros,utils.safeAdd,utils.safeNormalize]}])

    utils.showPlot()

if(__name__ == '__main__'):
    if(len(sys.argv)>1 and sys.argv[1]=='compare'):
        runComparison()
    else:
        runTestCases()
