#coding: utf-8
import os
import math
import utils
import pylab as pl
import numpy as np
import multiprocessing as mp
import otero_precipitation as op
from scipy.optimize import minimize
import datetime
import scipy.stats as stats

def calculateMetrics(time_range,RES,real_values):
    lwE=np.array([RES[(np.abs(time_range-t)).argmin(),op.EGG]-RES[(np.abs(time_range-(t-7))).argmin(),op.EGG] for t in time_range])
    lwE[lwE<0]=0.#replace negatives with zeros
    lwE=lwE/np.max(lwE)#normalize
    error=sum([(real_values[t]-lwE[t] )**2  for t in range(0,len(real_values)) if real_values[t]] )
    rho,p_value=stats.spearmanr(real_values[real_values!=[None]],lwE[real_values!=[None]])
    return lwE,error,rho,p_value

def error(x,real_values):
    op.vBS_od=x[0:op.n]#this is
    op.vBS_id=x[op.n:op.n+op.m]#a bad
    op.ws_s=x[op.n+op.m]#practice
    time_range,INPUT,RES = op.solveEquations()
    lwE,error,rho,p_value=calculateMetrics(time_range,RES,real_values)
    print('pid %i :'%os.getpid()),
    print(op.vBS_od),
    print(op.vBS_id),
    print(op.ws_s),
    print('Error: %s'%(error))
    return error

def getOptimalParameters(vNormalized_ovitrap_eggs):
    #initial value
    x0=np.append( np.append(op.vBS_od,op.vBS_id),np.array(op.ws_s))
    #Σ vBS_od[i] + vBS_id[i] = 1
    constraints = ({'type': 'eq', 'fun': lambda x:  1 - sum(x[0:op.n+op.m])})
    #0<=x<=1,0<=ws_s<=1.
    bounds=tuple((0,1) for x in x0 )#tuple((0,1) for x in range(0,len(x0)-1) ) + tuple((0,1.0) for x in range(0,1))

    opt=minimize(error,x0,(vNormalized_ovitrap_eggs),method='SLSQP',bounds=bounds,constraints=constraints,options={'eps': 1e-02, 'ftol': 1e-01})
    #print(res)
    return opt

def populate(time_range,ovitrap_eggs):
    result=np.array([None]*len(time_range))
    for d,eggs in enumerate(ovitrap_eggs):
        result[(np.abs(time_range-d)).argmin()]=eggs
    return result


def print_parameters(vBS_od,vBS_id,ws_s):
    print('ws_s=%f'%(ws_s))
    print("BS_o     * vo=%f * "%(round(sum(vBS_od),2)) + str( np.round( (1./sum(vBS_od))*vBS_od,2) ) )
    print("(1-BS_o) * vi=%f * "%(round(sum(vBS_id),2)) + str( np.round( (1./sum(vBS_id))*vBS_id,2) ) )
    return ''


if(__name__ == '__main__'):
    vNormalized_ovitrap_eggs=[]
    for ovitrap in range(1,30):
        ovitrap_eggs=utils.getOvitrapEggsFromCsv('data/private/Datos sensores de oviposicion.NO.csv',op.start_date,op.end_date,ovitrap)
        normalized_ovitrap_eggs=[e/max(ovitrap_eggs) if e else None for e in ovitrap_eggs]
        vNormalized_ovitrap_eggs.extend([populate(op.getTimeRange(),normalized_ovitrap_eggs)])

    print('initial Conditions:\n'),print_parameters(op.vBS_od,op.vBS_id,op.ws_s)
    p = mp.Pool(mp.cpu_count() -1)
    vOpt=p.map(getOptimalParameters, vNormalized_ovitrap_eggs)

    for i,opt in enumerate(vOpt):
        x=np.array(opt.x)
        op.vBS_od,op.vBS_id,op.ws_s=x[0:op.n],x[op.n:op.n+op.m],x[op.n+op.m]
        time_range,INPUT,RES = op.solveEquations()
        j=i%5
        if(j==0): pl.figure()
        pl.subplot(511 +j)
        ovitrap_eggs=utils.getOvitrapEggsFromCsv('data/private/Datos sensores de oviposicion.NO.csv',op.start_date,op.end_date,i+1)
        date_range=[datetime.timedelta(days=d)+datetime.datetime.combine(op.start_date,datetime.time()) for d in time_range]
        normalized_ovitrap_eggs=vNormalized_ovitrap_eggs[i]
        pl.plot(date_range, normalized_ovitrap_eggs, '^', label='Ovitrap %s'%(i+1),clip_on=False, zorder=100,markersize=8)

        lwE,error,rho,p_value=calculateMetrics(time_range,RES,normalized_ovitrap_eggs)
        pl.plot(date_range,lwE, '-k', label='Eggs normalized. Error: %s Rho:%s p-value:%s'%(error,rho,p_value))
        pl.legend(loc=0,prop={'size':10})

        print('Ovitrap %s:'%(i+1))
        print_parameters(op.vBS_od,op.vBS_id,op.ws_s)
        print(vOpt[i])
        print('-'*200)

    pl.xlabel('Time(in days starting in July)')
    pl.ylabel('')
    raw_input('Press enter')
    pl.show()
