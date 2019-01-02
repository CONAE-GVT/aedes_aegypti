#coding: utf-8
import os
import utils
import datetime
import numpy as np
import scipy.stats as stats
import multiprocessing as mp
from config import Configuration
from otero_precipitation import Model
from scipy.optimize import minimize,differential_evolution

MINIMIZE_METHOD='SLSQP'
OVITRAP_FILENAME='data/private/ovitrampas_2017-2018.csv'

def safeAdd(values):
    if(len(values.shape)!=2):
        return values
    else:
        return np.sum(values,axis=1)

def getConfiguration(x=None,n=None):
    if(x is None and n is None): return Configuration('resources/otero_precipitation.cfg')#no arguments passed, return default configuration

    x[0:n]/=x[0:n].sum()#sometimes the constraints fails
    start_date,end_date=utils.getStartEndDates(OVITRAP_FILENAME)
    configuration=Configuration('resources/otero_precipitation.cfg',
        {
        'breeding_site':{
            'amount':x[n],
            'distribution':x[0:n]
            },
        'simulation':{
            'start_date':start_date,
            'end_date':end_date
            }
        })
    return configuration

def calculateMetrics(time_range,mEggs,real_values):
    cEggs=safeAdd(mEggs)
    lwE=np.array([cEggs[(np.abs(time_range-t)).argmin()]-cEggs[(np.abs(time_range-(t-7))).argmin()] for t in time_range])
    lwE[lwE<0]=0.#replace negatives with zeros
    error=sum([(real_values[t]-lwE[t] )**2  for t in range(0,len(real_values)) if real_values[t]] )
    rho,p_value=stats.pearsonr(real_values[real_values!=[None]],lwE[real_values!=[None]])
    return lwE,error,rho,p_value

def populate(time_range,ovitrap_eggs):
    result=np.array([None]*len(time_range))
    for d,eggs in enumerate(ovitrap_eggs):
        result[(np.abs(time_range-d)).argmin()]=eggs
    return result

def error(x,ovitrap_eggs_i_with_id):
    #return np.dot(x,x)
    ovitrap_id,ovitrap_eggs_i=ovitrap_eggs_i_with_id
    model=Model(getConfiguration(x,len(x)-1))#vBS_d,BS_a#TODO:not agnostic
    ovitrap_eggs_i=populate(model.time_range,ovitrap_eggs_i)
    time_range,INPUT,Y=model.solveEquations()
    lwE,error,rho,p_value=calculateMetrics(time_range,Y[:,model.parameters.EGG],ovitrap_eggs_i)
    print('ovitrap:%s pid:%i d:%s BS_a:%s Error:%s rho:%s p value:%s'%(ovitrap_id,os.getpid(),model.parameters.vBS_d.tolist(),model.parameters.BS_a,error,rho,p_value))
    return error

def getOptimalParameters(ovitrap_eggs_i_with_id):
    configuration=getConfiguration()
    #initial value
    x0=np.append( configuration.getArray('breeding_site','distribution'), configuration.getFloat('breeding_site','amount'))
    #Σ vBS_d[i] = 1
    constraints = ({'type': 'eq', 'fun': lambda x:  1 - sum(x[0:-1])})#TODO:not agnostic
    #0<=x<=1,0<=ws_s<=1.
    bounds=tuple((1e-8,1) for x in x0[0:-1] )+ ((0,200),)#TODO:not agnostic

    if(MINIMIZE_METHOD=='SLSQP'):
        opt=minimize(error,x0,ovitrap_eggs_i_with_id,method='SLSQP',bounds=bounds,constraints=constraints,options={'eps': 1e-02, 'ftol': 1e-01})
    else:
        opt=differential_evolution(error,bounds,args=(ovitrap_eggs_i_with_id,))

    return opt

if(__name__ == '__main__'):
    start_date,end_date=utils.getStartEndDates(OVITRAP_FILENAME)
    ovitrap_eggs=[[ovitrap_id,utils.getOvitrapEggsFromCsv(OVITRAP_FILENAME,start_date,end_date,ovitrap_id)] for ovitrap_id in range(2,151)]

    print('Starting...')
    vOpt=mp.Pool(mp.cpu_count()-2).map(getOptimalParameters, ovitrap_eggs)
    print(vOpt)
