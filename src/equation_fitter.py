#coding: utf-8
import os
import math
import utils
import datetime
import pylab as pl
import numpy as np
import scipy.stats as stats
import multiprocessing as mp
import otero_precipitation as op
import performance_equations as pe
from scipy.optimize import minimize,differential_evolution

MINIMIZE_METHOD='differential_evolution'

def calculateMetrics(time_range,RES,real_values):
    lwE=np.array([RES[(np.abs(time_range-t)).argmin(),op.EGG]-RES[(np.abs(time_range-(t-7))).argmin(),op.EGG] for t in time_range])
    lwE[lwE<0]=0.#replace negatives with zeros
    lwE=lwE/np.max(lwE)#normalize
    error=sum([(real_values[t]-lwE[t] )**2  for t in range(0,len(real_values)) if real_values[t]] )
    rho,p_value=stats.spearmanr(real_values[real_values!=[None]],lwE[real_values!=[None]])
    return lwE,error,rho,p_value

def error(x,real_values,ovitrap=None):
    if(not ovitrap):
        real_values,ovitrap=real_values#weird differential_evolution bug...

    l=np.sum(x[0:op.n+op.m])
    if(l<1e-5): return 500.
    x[0:op.n+op.m]/=l#constraint: #Σ vBS_od[i] + vBS_id[i] = 1

    op.vBS_od=x[0:op.n]#this is
    op.vBS_id=x[op.n:op.n+op.m]#a bad
    op.ws_s=x[op.n+op.m]#practice
    MAX_BS_A=4550.
    op.BS_a=MAX_BS_A*x[op.n+op.m+1]+50#practice
    time_range,INPUT,RES = pe.solvePopulationEquations(initial_condition = [100.0, 0.0,0.0,0.0,0.0])
    lwE,error,rho,p_value=calculateMetrics(time_range,RES,real_values)
    ovitrap_data_len=float(len(real_values[real_values!=[None]]))
    if(ovitrap_data_len==0): ovitrap_data_len=1.
    print('ovitrap:%s pid:%i od:%s id:%s ws_s:%s BS_a:%s Error:%s len:%s Error/len: %s'%(ovitrap,os.getpid(),op.vBS_od,op.vBS_id,op.ws_s,op.BS_a,error,ovitrap_data_len,error/ovitrap_data_len))
    return error

def earlyOut(xk,convergence):
    print('*'*20)
    if(convergence>0.5): return True

def getOptimalParameters(args):
    #initial value
    x0=np.append( np.append(op.vBS_od,op.vBS_id),np.array([op.ws_s,op.BS_a]))
    #Σ vBS_od[i] + vBS_id[i] = 1
    constraints = ({'type': 'eq', 'fun': lambda x:  1 - sum(x[0:op.n+op.m])})
    #0<=x<=1,0<=ws_s<=1.
    bounds=tuple((0,1) for x in x0 )#tuple((0,1) for x in range(0,len(x0)-1) ) + tuple((0,1.0) for x in range(0,1))

    opt=None
    if(MINIMIZE_METHOD=='SLSQP'):
        opt=minimize(error,x0,args,method='SLSQP',bounds=bounds,constraints=constraints,options={'eps': 1e-02, 'ftol': 1e-01})
    else:
        opt=differential_evolution(error,bounds,args=args,callback=earlyOut,maxiter=1, popsize=1)

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
    time_range,W0,W=pe.solveWaterEquations(initial_condition=[0. for i in op.vBS_od])
    vNormalized_ovitrap_eggs=[]
    for ovitrap in range(1,30):
        ovitrap_eggs=utils.getOvitrapEggsFromCsv('data/private/Datos sensores de oviposicion.NO.csv',op.start_date,op.end_date,ovitrap)
        ovitrap_eggs=np.array(ovitrap_eggs)
        normalized_ovitrap_eggs=[e/ovitrap_eggs[np.not_equal(ovitrap_eggs,None)].max() if e else None for e in ovitrap_eggs]
        vNormalized_ovitrap_eggs.extend([[populate(op.getTimeRange(),ovitrap_eggs),ovitrap]])

    print('initial Conditions:\n'),print_parameters(op.vBS_od,op.vBS_id,op.ws_s)
    p = mp.Pool(mp.cpu_count() -2)
    vOpt=p.map(getOptimalParameters, vNormalized_ovitrap_eggs)

    #
    weight= lambda z: 1./z
    total_weight=np.sum([weight(opt.fun+1e-1) for opt in vOpt])
    weight_mean_x=np.array([weight(opt.fun+1e-1)*opt.x/total_weight for opt in vOpt]).sum(axis=0)
    print(weight_mean_x)
    for idx,opt in enumerate(vOpt):
        x=np.array(opt.x)
        op.vBS_od,op.vBS_id,op.ws_s=x[0:op.n],x[op.n:op.n+op.m],x[op.n+op.m]
        time_range,INPUT,RES = op.solveEquations(initial_condition = [100.0, 0.0,0.0,0.0,0.0]+[0 for i in range(0,op.n)])
        j=idx%5
        if(j==0): pl.figure()
        pl.subplot(511 +j)
        ovitrap_eggs=utils.getOvitrapEggsFromCsv('data/private/Datos sensores de oviposicion.NO.csv',op.start_date,op.end_date,idx+1)
        date_range=[datetime.timedelta(days=d)+datetime.datetime.combine(op.start_date,datetime.time()) for d in time_range]
        normalized_ovitrap_eggs=vNormalized_ovitrap_eggs[idx]
        pl.plot(date_range, normalized_ovitrap_eggs, '^', label='Ovitrap %s'%(idx+1),clip_on=False, zorder=100,markersize=8)

        lwE,error,rho,p_value=calculateMetrics(time_range,RES,normalized_ovitrap_eggs)
        pl.plot(date_range,lwE, '-k', label='Eggs normalized. Error: %s Rho:%s p-value:%s'%(error,rho,p_value))
        pl.legend(loc=0,prop={'size':10})

        print('Ovitrap %s:'%(idx+1))
        print_parameters(op.vBS_od,op.vBS_id,op.ws_s)
        print(vOpt[idx])
        print('-'*200)


    pl.xlabel('Time(in days starting in July)')
    pl.ylabel('')
    raw_input('Press enter')
    pl.show()
