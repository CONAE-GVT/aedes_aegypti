#coding: utf-8
import os
import math
import utils
import pylab as pl
import numpy as np
import multiprocessing as mp
import otero_precipitation as op
from scipy.optimize import minimize

def nL_error(x,real_values):#Normalized Larvaes error#TODO: maybe it's best to use an interpolator polynomial of the real values
    op.vBS_od=x[0:op.n]#this is
    op.vBS_id=x[op.n:op.n+op.m]#a bad
    op.ws_s=x[op.n+op.m]#practice
    time_range,INPUT,RES = op.solveEquations()
    nL=(1.0/max(RES[:,op.LARVAE])) * np.array(RES[:,op.LARVAE])
    error=sum([(real_values[t]-nL[t] )**2  for t in range(0,len(real_values)) if real_values[t]] )
    print('pid %i :'%os.getpid()),
    print(op.vBS_od),
    print(op.vBS_id),
    print(op.ws_s),
    print('Error: %s'%(error))
    return error

def getOptimalParameters(normalized_ris):
    #initial value
    x0=np.append( np.append(op.vBS_od,op.vBS_id),np.array(op.ws_s))
    #Σ vBS_od[i] + vBS_id[i] = 1
    constraints = ({'type': 'eq', 'fun': lambda x:  1 - sum(x[0:op.n+op.m])})
    #0<=x<=1,0<=ws_s<=1.
    bounds=tuple((0,1) for x in x0 )#tuple((0,1) for x in range(0,len(x0)-1) ) + tuple((0,1.0) for x in range(0,1))

    opt=minimize(nL_error,x0,(normalized_ris),method='SLSQP',bounds=bounds,constraints=constraints)
    #print(res)
    return opt

def populate(indices,values):#return an array containing the i-eth value in the i-eth  indices
    if(not len(indices)==len(values)):
        raise ValueError('indices and values must have the same length')
    indices=np.array(indices)
    values=np.array(values)
    return [values[ np.where(indices==i)[0][0] ] if len(np.where(indices==i)[0])==1 else None for i in range(0,max(indices)+1) ]

def print_parameters(vBS_od,vBS_id,ws_s):
    print('ws_s=%f'%(ws_s))
    print("BS_o     * vo=%f * "%(round(sum(vBS_od),2)) + str( np.round( (1./sum(vBS_od))*vBS_od,2) ) )
    print("(1-BS_o) * vi=%f * "%(round(sum(vBS_id),2)) + str( np.round( (1./sum(vBS_id))*vBS_id,2) ) )
    return ''


if(__name__ == '__main__'):
    vHi_days=[]
    vZones=[]
    vNormalized_his=[]
    for zone in op.location['zones']:
        hi_days, his= utils.getIndexesForPlot(op.AEDIC_INDICES_FILENAME,op.start_date,0,2,zone,1)
        normalized_his=np.multiply(1.0/max(his),his)
        vZones.extend([zone])#just for plotting
        vHi_days.extend([hi_days])#just for plotting
        vNormalized_his.extend([populate(hi_days,normalized_his)])
    print('initial Conditions:\n'),print_parameters(op.vBS_od,op.vBS_id,op.ws_s)
    p = mp.Pool(mp.cpu_count() -1)
    vOpt=p.map(getOptimalParameters, vNormalized_his)

    for i,opt in enumerate(vOpt):
        x=np.array(opt.x)
        op.vBS_od,op.vBS_id,op.ws_s=x[0:op.n],x[op.n:op.n+op.m],x[op.n+op.m]
        time_range,INPUT,RES = op.solveEquations()
        pl.subplot(511 +i)
        pl.plot(time_range, np.multiply(1.0/max(RES[:,op.LARVAE]), RES[:,op.LARVAE]), '-', label=vZones[i] +' Larvaes normalized')
        pl.plot(range(0,len(vNormalized_his[i])),vNormalized_his[i], '^', label=vZones[i] +' Larvaes  normalized',clip_on=False, zorder=100,markersize=8)
        pl.legend(loc=0,prop={'size':10})

        print(vZones[i])
        print_parameters(op.vBS_od,op.vBS_id,op.ws_s)
        print(vOpt[i])
        print('-'*200)

    pl.xlabel('Time(in days starting in July)')
    pl.ylabel('')
    raw_input('Press enter')
    pl.show()
