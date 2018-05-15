#coding: utf-8
import os
import math
import utils
import datetime
import pylab as pl
import numpy as np
import scipy.stats as stats
import multiprocessing as mp
from _equations import diff_eqs
from config import Configuration
from otero_precipitation import Model
#import performance_equations as pe
from scipy.optimize import minimize,differential_evolution

MINIMIZE_METHOD='differential_evolution'

def safeAdd(values):
    if(len(values.shape)!=2):
        return values
    else:
        return np.sum(values,axis=1)

def calculateMetrics(time_range,mEggs,real_values):
    cEggs=safeAdd(mEggs)
    lwE=np.array([cEggs[(np.abs(time_range-t)).argmin()]-cEggs[(np.abs(time_range-(t-7))).argmin()] for t in time_range])
    lwE[lwE<0]=0.#replace negatives with zeros
    lwE=lwE/np.max(lwE)#normalize
    error=sum([(real_values[t]-lwE[t] )**2  for t in range(0,len(real_values)) if real_values[t]] )
    rho,p_value=stats.spearmanr(real_values[real_values!=[None]],lwE[real_values!=[None]])
    return lwE,error,rho,p_value

def getConfiguration(x,n,m):
    l=np.sum(x[0:n+m])
    if(l<1e-5): return 500.
    x[0:n+m]/=l#constraint: #Σ vBS_od[i] + vBS_id[i] = 1
    MAX_BS_A=4550.
    configuration=Configuration('resources/otero_precipitation.cfg',
        {'breeding_site':{
            'amount':MAX_BS_A*x[n+m]+50,
            'outside_distribution':x[0:n],
            'inside_distribution':x[n:n+m]
            }
        })
    return configuration

def error(x,model,real_values=None,ovitrap=None):
    #return np.dot(x,x)
    if(real_values is None):
        model,real_values,ovitrap=model#weird differential_evolution bug...

    model=Model(getConfiguration(x,model.parameters.n,model.parameters.m))
    #time_range,INPUT,RES = pe.solvePopulationEquations(initial_condition = [100.0, 0.0,0.0,0.0,0.0])
    time_range,INPUT,RES=model.solveEquations(equations=diff_eqs,method='rk')
    lwE,error,rho,p_value=calculateMetrics(time_range,RES[:,model.parameters.EGG],real_values)
    ovitrap_data_len=float(len(real_values[real_values!=[None]]))
    if(ovitrap_data_len==0): ovitrap_data_len=1.
    print('ovitrap:%s pid:%i od:%s id:%s BS_a:%s Error:%s len:%s Error/len: %s'%(ovitrap,os.getpid(),model.parameters.vBS_od,model.parameters.vBS_id,model.parameters.BS_a,error,ovitrap_data_len,error/ovitrap_data_len))
    return error

def getOptimalParameters(args):
    model=Model()
    #initial value
    x0=np.append( np.append(model.parameters.vBS_od,model.parameters.vBS_id),np.array([model.parameters.BS_a]))
    #Σ vBS_od[i] + vBS_id[i] = 1
    constraints = ({'type': 'eq', 'fun': lambda x:  1 - sum(x[0:model.parameters.n+model.parameters.m])})
    #0<=x<=1,0<=ws_s<=1.
    bounds=tuple((0,1) for x in x0 )#tuple((0,1) for x in range(0,len(x0)-1) ) + tuple((0,1.0) for x in range(0,1))

    opt=None
    args=[model]+args
    if(MINIMIZE_METHOD=='SLSQP'):
        opt=minimize(error,x0,args,method='SLSQP',bounds=bounds,constraints=constraints,options={'eps': 1e-02, 'ftol': 1e-01})
    else:
        opt=differential_evolution(error,bounds,args=args)

    return opt

def populate(time_range,ovitrap_eggs):
    result=np.array([None]*len(time_range))
    for d,eggs in enumerate(ovitrap_eggs):
        result[(np.abs(time_range-d)).argmin()]=eggs
    return result


def print_parameters(vBS_od,vBS_id):
    print("BS_o     * vo=%f * "%(round(sum(vBS_od),2)) + str( np.round( (1./sum(vBS_od))*vBS_od,2) ) )
    print("(1-BS_o) * vi=%f * "%(round(sum(vBS_id),2)) + str( np.round( (1./sum(vBS_id))*vBS_id,2) ) )
    return ''

if(__name__ == '__main__'):
    #time_range,W0,W=pe.solveWaterEquations(initial_condition=[0. for i in op.vBS_od])
    vNormalized_ovitrap_eggs=[]
    for ovitrap in range(1,30):
        model=Model()#just for dates and time_range
        ovitrap_eggs=utils.getOvitrapEggsFromCsv('data/private/Datos sensores de oviposicion.NO.csv',model.start_date,model.end_date,ovitrap)
        ovitrap_eggs=np.array(ovitrap_eggs)
        normalized_ovitrap_eggs=[e/ovitrap_eggs[np.not_equal(ovitrap_eggs,None)].max() if e else None for e in ovitrap_eggs]
        vNormalized_ovitrap_eggs.extend([[populate(model.time_range,ovitrap_eggs),ovitrap]])

    print('Starting...')
    p = mp.Pool(mp.cpu_count() -2)
    vOpt=p.map(getOptimalParameters, vNormalized_ovitrap_eggs)

    #
    weight= lambda z: 1./z
    total_weight=np.sum([weight(opt.fun+1e-1) for opt in vOpt])
    weight_mean_x=np.array([weight(opt.fun+1e-1)*opt.x/total_weight for opt in vOpt]).sum(axis=0)
    print(weight_mean_x)
    for idx,opt in enumerate(vOpt):
        x=np.array(opt.x)
        model=Model()
        n,m=model.parameters.n,model.parameters.m
        model.parameters.vBS_od,model.parameters.vBS_id=x[0:n],x[n:n+m]
        time_range,INPUT,RES = model.solveEquations(equations=diff_eqs,method='rk')
        filename=model.save()
        new_filename=filename.replace('.csv','_ovi%s.csv'%(idx+1))
        os.rename(filename,new_filename)
        os.rename(filename.replace('.csv','.cfg'),new_filename.replace('.csv','.cfg'))
        print('Ovitrap %s:'%(idx+1))
        print_parameters(model.parameters.vBS_od,model.parameters.vBS_id)
        print(vOpt[idx])
        print('-'*200)
