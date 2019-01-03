#coding: utf-8
import os
import utils
import datetime
import pylab as pl
import numpy as np
import scipy.stats as stats
import multiprocessing as mp
from config import Configuration
from otero_precipitation import Model
from scipy.optimize import minimize,differential_evolution

MINIMIZE_METHOD='SLSQP'
OVITRAP_FILENAME='data/private/ovitrampas_2017-2018.csv'

class perOvitrap:
    def __init__(self,model):
        self.BS_a=model.parameters.BS_a
        self.vBS_d=model.parameters.vBS_d

    def __call__(self,lwE):
        if (lwE.ndim!=2):
            return  lwE#if ndims!=2 is not the lwE, but the real ovitraps#not a clean way to check this.
        else:
            return lwE[:,0]/(self.BS_a*self.vBS_d[0])

def getConfiguration(x=None,n=None):
    if(x is None and n is None): return Configuration('resources/otero_precipitation.cfg')#no arguments passed, return default configuration

    start_date,end_date=utils.getStartEndDates(OVITRAP_FILENAME)
    configuration=Configuration('resources/otero_precipitation.cfg',
        {
        'breeding_site':{
            'manually_filled':x[0:n]
            },
        'simulation':{
            'start_date':start_date,
            'end_date':end_date
            },
        'biology':{
            'alpha0':x[n:]
            }
        })
    return configuration

def calculateMetrics(time_range,model,ovitrap_eggs_i):
    mEggs=model.Y[:,model.parameters.EGG]
    lwE=np.array([mEggs[(np.abs(time_range-t)).argmin()]-mEggs[(np.abs(time_range-(t-7))).argmin()] for t in time_range])#last week eggs
    lwE[lwE<0]=0.#replace negatives with zeros
    lwE=perOvitrap(model)(lwE)
    error=sum([(ovitrap_eggs_i[t]-lwE[t] )**2  for t in range(0,len(ovitrap_eggs_i)) if ovitrap_eggs_i[t]] )
    rho,p_value=stats.pearsonr(ovitrap_eggs_i[ovitrap_eggs_i!=[None]],lwE[ovitrap_eggs_i!=[None]])
    return lwE,error,rho,p_value

def populate(time_range,ovitrap_eggs):
    result=np.array([None]*len(time_range))
    for d,eggs in enumerate(ovitrap_eggs):
        result[(np.abs(time_range-d)).argmin()]=eggs
    return result

def error(x,ovitrap_eggs_i_with_id):
    #return np.dot(x,x)
    ovitrap_id,ovitrap_eggs_i=ovitrap_eggs_i_with_id
    model=Model(getConfiguration(x,int(len(x)/2)))#vBS_d,BS_a#TODO:not agnostic
    ovitrap_eggs_i=populate(model.time_range,ovitrap_eggs_i)
    time_range,INPUT,Y=model.solveEquations(method='rk')
    lwE,error,rho,p_value=calculateMetrics(time_range,model,ovitrap_eggs_i)
    title='ovitrap:%s\nmf:%s\n' r'$\alpha_0$:%s' '\n' r'Error:%s $\rho$:%s p value:%s'%(ovitrap_id,model.parameters.vBS_mf.tolist(),model.parameters.vAlpha0.tolist(),error,rho,p_value)
    utils.plot(model,subplots=[{'lwE':'','f':[utils.replaceNegativesWithZeros,perOvitrap(model)],'O':[int(ovitrap_id)]}],title=title,figure=False)
    return error

def getOptimalParameters(ovitrap_eggs_i_with_id):
    pl.figure()#hack
    configuration=getConfiguration()
    vmf=configuration.getArray('breeding_site','manually_filled')
    vAlpha0=configuration.getArray('biology','alpha0')
    #initial value
    x0=np.append( vmf, vAlpha0)
    #Σ vBS_d[i] = 1
    constraints = ()#({'type': 'eq', 'fun': lambda x:  1 - sum(x[0:-1])})#TODO:not agnostic
    #0<=x<=1,0<=ws_s<=1.
    bounds=tuple((1e-8,1) for x in vmf )+ tuple((1e-8,1) for x in vAlpha0)#TODO:not agnostic

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
