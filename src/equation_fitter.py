#coding: utf-8
import os
import sys
import utils
import datetime
import numpy as np
import scipy.stats as stats
import multiprocessing as mp
from config import Configuration
from otero_precipitation_wrapper_wrapper import Model
from scipy.optimize import minimize,differential_evolution

MINIMIZE_METHOD='differential_evolution'
OVITRAP_FILENAME='data/private/ovitrampas_2017-2019.mean.csv'

def getConfiguration(x=None,n=None):
    CONFIGURATION_FILENAME='resources/1c.cfg'
    if(x is None and n is None): return Configuration(CONFIGURATION_FILENAME)#no arguments passed, return default configuration

    start_date,end_date=utils.getStartEndDates(OVITRAP_FILENAME)
    configuration=Configuration(CONFIGURATION_FILENAME,
        {
        'breeding_site':{
            'height':x[0:n],
            'manually_filled':x[n:2*n],
            'bare':x[2*n:3*n],
            'evaporation_factor':x[3*n:]
            },
        'simulation':{
            #'start_date':start_date,#TODO:start date should be fixed, not the same as ovitrap start date
            'end_date':end_date
            }
        })
    return configuration

def calculateMetrics(time_range,model,ovitrap_eggs_i):
    BS_a,vBS_d,m,n,OVIPOSITION=model.parameters.BS_a,model.parameters.vBS_d,model.parameters.m,model.parameters.n,model.parameters.OVIPOSITION
    Y=model.Y
    indexOf=lambda t: (np.abs(time_range-t)).argmin()
    lwO=np.array([ (Y[indexOf(t),OVIPOSITION]-Y[indexOf(t-7),OVIPOSITION]).reshape(m,n).sum(axis=0) for t in time_range])/(BS_a*vBS_d)
    lwO=lwO[:,0]#if multiple container w assume the first one is ovitrap and the rest are wild containers
    d=utils.rmse(ovitrap_eggs_i[ovitrap_eggs_i!=[None]], lwO[ovitrap_eggs_i!=[None]])
    return d

def populate(time_range,start_date,ovitrap_eggs):
    result=np.array([None]*len(time_range))
    for date in ovitrap_eggs:
        d=(date-start_date).days
        eggs=utils.noneMean(ovitrap_eggs[date])
        result[(np.abs(time_range-d)).argmin()]=eggs
    return result

def error(x,ovitrap_eggs_i_with_id):
    #return np.dot(x,x)
    ovitrap_id,ovitrap_eggs_i=ovitrap_eggs_i_with_id
    model=Model(getConfiguration(x,int(len(x)/4)))#TODO:not agnostic
    ovitrap_eggs_i=populate(model.time_range,model.start_date,ovitrap_eggs_i)
    time_range,Y=model.solveEquations()
    error=calculateMetrics(time_range,model,ovitrap_eggs_i)
    #title='ovitrap:%s\nmf:%s\n' r'$\alpha_0$:%s' '\n' r'Error:%s $\rho$:%s p value:%s'%(ovitrap_id,model.parameters.vBS_mf.tolist(),model.parameters.vAlpha0.tolist(),error,rho,p_value)
    #print(title)
    return error

def getOptimalParameters(ovitrap_eggs_i_with_id):
    configuration=getConfiguration()
    vBS_h=configuration.getArray('breeding_site','height')
    vmf=configuration.getArray('breeding_site','manually_filled')
    vBS_b=configuration.getArray('breeding_site','bare')
    vBS_ef=configuration.getArray('breeding_site','evaporation_factor')
    #initial value
    x0=np.concatenate( (vBS_h, vmf,vBS_b,vBS_ef) )
    #Σ vBS_d[i] = 1
    constraints = ()#({'type': 'eq', 'fun': lambda x:  1 - sum(x[0:-1])})#TODO:not agnostic
    #0<=x<=1,0<=ws_s<=1.
    bounds=tuple((2,20) for x in vBS_h) + tuple((0,2) for x in vmf ) + tuple((0,1) for x in vBS_b) + tuple((0,2) for x in vBS_ef)#TODO:not agnostic

    if(MINIMIZE_METHOD=='SLSQP'):
        opt=minimize(error,x0,ovitrap_eggs_i_with_id,method='SLSQP',bounds=bounds,constraints=constraints,options={'eps': 1e-02, 'ftol': 1e-01})
    else:
        opt=differential_evolution(error,bounds,args=(ovitrap_eggs_i_with_id,))#, workers=mp.cpu_count()-2 for scipy>1.3

    open('out/equation_fitter/_ovi%s.txt'%ovitrap_eggs_i_with_id[0],'w').write(str(opt)+'\n\n')#ovitrap_eggs_i_with_id[0] is the ovitrap_id

    return opt

if(__name__ == '__main__'):
    start_date,end_date=utils.getStartEndDates(OVITRAP_FILENAME)
    if(len(sys.argv)>1):
        ovitrap_id=int(sys.argv[1])
        ovitrap_eggs_i_with_id=[ovitrap_id,utils.getOvitrapEggsFromCsv2(OVITRAP_FILENAME,start_date,end_date,ovitrap_id)]
        vOpt=getOptimalParameters(ovitrap_eggs_i_with_id)
    else:
        ovitrap_eggs=[[ovitrap_id,utils.getOvitrapEggsFromCsv2(OVITRAP_FILENAME,start_date,end_date,ovitrap_id)] for ovitrap_id in range(1,152)]

        print('Starting...')
        vOpt=mp.Pool(mp.cpu_count()-2).map(getOptimalParameters, ovitrap_eggs)
        # vOpt=[]
        # for i,ovitrap_eggs_i_with_id in enumerate(ovitrap_eggs):
        #     vOpt+=getOptimalParameters(ovitrap_eggs_i_with_id)
        #     print('\r %s'%(round(float(i)/len(ovitrap_eggs)*100,2) ), end='')
        print(vOpt)
