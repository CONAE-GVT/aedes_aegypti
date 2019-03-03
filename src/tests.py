#coding: utf-8
import os
import re
import sys
import utils
import datetime
import numpy as np
from config import Configuration
from otero_precipitation import Model
from equations import diff_eqs,vR_D
import equations
from spatial_equations import diff_eqs as spatial_diff_eqs
import pylab as pl

def testModel(configuration, subplots=[],plot_start_date=None,title='',figure=True,color=None):
    model=Model(configuration)
    time_range,INPUT,RES=model.solveEquations(equations=utils.MetricsEquations(model,diff_eqs),method='rk' )
    utils.plot(model,subplots,plot_start_date,title,figure,color)
    for warning in model.warnings:
        print('# WARNING: ' + warning)

def runOviShow(folder):
    config_filenames=[filename for filename in os.listdir(folder) if filename.endswith('.cfg')]
    for config_filename in config_filenames:
        for ovitrap_id in re.findall(r'.*_ovi([0-9]+)\.cfg',config_filename):
            testModel(Configuration(folder+'/'+config_filename),subplots=[ ['E','P','A1+A2',[utils.safeAdd,utils.normalize] ],{'lwE':'','O':[int(ovitrap_id)+1],'f':[utils.replaceNegativesWithZeros,utils.safeAdd,utils.safeNormalize]}])


def runShow(folder):
    filenames=[filename for filename in os.listdir(folder) if re.match(r'.*cordoba\.[0-9]+.*\.csv',filename)]
    filenames.sort()
    for i,filename in enumerate(filenames):
        configuration=Configuration('resources/otero_precipitation.cfg',
            {'simulation':{
                'start_date':datetime.date(2017,7,1),
                'end_date':datetime.date(2019,1,4),
            }
            })
        configuration.config_parser.set('location','name',filename.replace('.csv',''))
        p=i/len(filenames)
        color=p*np.array([1,0,0]) + (1-p)*np.array([0,1,0])
        testModel(configuration,subplots=[ ['A1+A2',[utils.safeAdd] ],['T'] ],plot_start_date=datetime.date(2018,12,1),title=': '+configuration.getString('location','name'),color=color.tolist(),figure=False)


def runProject():
    config=Configuration('resources/otero_precipitation.cfg',
        {'simulation':{
            'start_date':datetime.date(2017,7,1),
            'end_date':datetime.date.today()+datetime.timedelta(30),
        }
        })
    for location in utils.getLocations():
        config.config_parser.set('location','name',location+'.full')
        testModel(config,subplots=[['E','A1+A2',[utils.safeAdd,utils.normalize] ]],title=location)


def runSpatial():
    configuration=Configuration('resources/otero_precipitation.cfg')
    configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(30)))
    model=Model(configuration)
    #modify some parameters to make them spatial
    parameters=model.parameters
    n=parameters.n
    parameters.FLYER=3*n+2#in R
    parameters.diff=configuration.getFloat('biology','diffusion')#diffusion-like coefficient
    parameters.P,warning=utils.getPreferenceMatrix('data/public/goodness/'+configuration.getString('biology','goodness'),patch_size=100)
    model.warnings.append(warning)
    HEIGHT,WIDTH=parameters.P.shape[:2]
    parameters.initial_condition=np.append(parameters.initial_condition,[0])#append flyers
    parameters.initial_condition=(parameters.initial_condition*utils.getY0FactorMatrix(HEIGHT,WIDTH)[:,:,np.newaxis]).reshape(HEIGHT*WIDTH*(3*n+3))#TODO:find a better way of introducing initial conditions to spatial
    parameters.vBS_a=parameters.BS_a*np.ones((HEIGHT,WIDTH))#TODO:Estimate BS_a as in section 8.3 Otero 2008, "Estimation of the breeding site density"
    parameters.vBS_d=parameters.vBS_d*np.ones((HEIGHT,WIDTH,n))#ASSUMPTION: distribuition is constant on (x,y)
    parameters.vAlpha=(parameters.vAlpha0*np.ones((HEIGHT,WIDTH,n)) )/parameters.vBS_a[:,:,np.newaxis]
    theta=parameters.vBS_a/150.
    tdep=0.229#Average time for egg deposition Christophers(1960)
    parameters.ovr= np.where(parameters.vBS_a<=150, theta/tdep, 1/tdep)

    for warning in model.warnings:
        print('# WARNING: ' + warning)

    #solve the equations
    time_range,initial_condition,Y=model.solveEquations(equations=utils.ProgressEquations(model,spatial_diff_eqs),method='cuda_rk' )
    Y=Y.reshape(Y.shape[0],HEIGHT,WIDTH,3*n + 3)
    np.save('out/Y.npy',Y)
    #time_range,Y=model.time_range,np.load('out/Y.npy')#to debug video

    EGG,LARVAE,PUPAE,ADULT1,FLYER,ADULT2=parameters.EGG,parameters.LARVAE,parameters.PUPAE,parameters.ADULT1,parameters.FLYER,parameters.ADULT2
    stages={'E':EGG, 'A':[ADULT1,FLYER,ADULT2]}
    for key in stages:
        print('Creating animation for %s...'%key)
        matrix=np.sum(Y[:,:,:,stages[key]],axis=3)
        matrix=matrix/matrix.max()
        start_date=configuration.getDate('simulation','start_date')
        getTitle=lambda i: datetime.timedelta(days=time_range[i])+start_date
        utils.createAnimation('out/%s'%key,matrix,getTitle,time_range.max())# 1 day : 1 second

import equation_fitter
import scipy.stats as stats
from scipy import interpolate
def calculateMetrics(lwO_mean,lwO_std,ovitrap_eggs_i):
    lwO_u=np.sum(lwO_mean+lwO_std,axis=1)
    lwO_l=np.sum(lwO_mean-lwO_std,axis=1)
    lwO_2u=np.sum(lwO_mean+1.5*lwO_std,axis=1)
    lwO_2l=np.sum(lwO_mean- 1.5*lwO_std,axis=1)
    V=np.sum(lwO_mean,axis=1)

    V/=np.max(V)
    ovitrap_eggs_i/=np.nanmax(ovitrap_eggs_i)
    d=np.nansum((ovitrap_eggs_i-V)**2 )/ np.where(~np.isnan(ovitrap_eggs_i),1,0).sum()

    rho,p_value=stats.pearsonr(ovitrap_eggs_i[~np.isnan(ovitrap_eggs_i)], V[~np.isnan(ovitrap_eggs_i)])

    valid_ovi_idx=~np.isnan(ovitrap_eggs_i)
    count=np.sum( np.where(np.logical_and(lwO_l[valid_ovi_idx]<=ovitrap_eggs_i[valid_ovi_idx],ovitrap_eggs_i[valid_ovi_idx]<=lwO_u[valid_ovi_idx]),1,0)  )
    count2=np.sum(np.where(np.logical_and(lwO_2l[valid_ovi_idx]<=ovitrap_eggs_i[valid_ovi_idx],ovitrap_eggs_i[valid_ovi_idx]<=lwO_2u[valid_ovi_idx]),1,0)  )
    return d,count,count2,rho,p_value

def runCases(case):
    if(case==0):
        ovi_range=range(1,151)
        for h in [10]:
            for mf  in np.arange(0,10,0.1):
                configuration=Configuration('resources/2c.cfg')
                configuration.config_parser.set('location','name','cordoba.full')#TODO:fix data and
                configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(30)))# uncomment these two
                n=len(configuration.getArray('breeding_site','height'))
                configuration.config_parser.set('breeding_site','height',','.join([str(h)]*n))
                configuration.config_parser.set('breeding_site','manually_filled',','.join([str(mf)]+[str(0)]*(n-1)))
                model=Model(configuration)
                time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')

                errors=[[[1e5,0,0,-1,1],[1e5,0,0,-1,1]]]*151#just to fill the ovitrap 0 that do not exist in reality
                for ovitrap_id in ovi_range:
                    OVITRAP_FILENAME='data/private/ovitrampas_2017-2018.full.csv'
                    values=utils.getOvitrapEggsFromCsv2(OVITRAP_FILENAME,None,None,ovitrap_id)
                    ovitrap_days=values.keys()
                    dates=[model.start_date + datetime.timedelta(t) for t in time_range]

                    ovi_a=[values[date][0] if date in values else None for date in dates]#TODO:WARNING!this will repeat values if model granularity is not 1 value per day.
                    ovi_a=np.array(equation_fitter.populate(model.time_range,ovi_a))
                    ovi_a=np.array(ovi_a,dtype=np.float)#this change None for np.nan

                    ovi_b=[values[date][1] if (date in values and len(values[date])>1) else None for date in dates]#TODO:WARNING!this will repeat values if model granularity is not 1 value per day.Also this is different than test 13
                    ##ovi_b=[values[date][1] if len(values[date])>1 else values[date][0] for date in ovitrap_days]
                    ovi_b=np.array(equation_fitter.populate(model.time_range,ovi_b))
                    ovi_b=np.array(ovi_b,dtype=np.float)#this change None for np.nan

                    indexOf=lambda t: (np.abs(time_range-t)).argmin()
                    OVIPOSITION=model.parameters.OVIPOSITION
                    BS_a=model.parameters.BS_a
                    O=Y[:,OVIPOSITION]
                    lwO=np.array([Y[indexOf(t),OVIPOSITION]-Y[indexOf(t-7),OVIPOSITION] for t in time_range])/BS_a
                    lwO_mean=np.array([lwO[indexOf(t-7):indexOf(t+7)].mean(axis=0) for t in time_range])
                    lwO_std =np.array([lwO[indexOf(t-7):indexOf(t+7)].std(axis=0) for t in time_range])



                    square_a,count_a,count2_a,rho_a,p_value_a=calculateMetrics(lwO_mean,lwO_std,ovi_a)
                    square_b,count_b,count2_b,rho_b,p_value_b=calculateMetrics(lwO_mean,lwO_std,ovi_b)
                    errors[ovitrap_id]=[[square_a,count_a,count2_a,rho_a,p_value_a],[square_b,count_b,count2_b,rho_b,p_value_b]]

                for i,ovi_type in enumerate(['a','b']):
                    error=np.array(errors)
                    error=error[:,i,:]
                    square,count,count2,rho=error[:,0],error[:,1],error[:,2],error[:,3]
                    sort0,sort1,sort2,sort3=np.argsort(square),np.argsort(count),np.argsort(count2),np.argsort(rho)
                    print('ovi:%s mf:%scm. h: %scm.---->(best) \t square: id:%s,score:%s \t count:id:%s,score:%s \t count2:id:%s,score:%s \t rho:id:%s,score:%s'%(ovi_type,mf,h,sort0[0],square[sort0[0]], sort1[-1],count[sort1[-1]], sort2[-1],count2[sort2[-1]],sort3[-1],rho[sort3[-1]] ) )
        pl.show()

    if(case==1):
        for mf  in [3,3.1,3.2,3.3]:
            h=10.
            configuration=Configuration('resources/2c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
            n=len(configuration.getArray('breeding_site','height'))
            configuration.config_parser.set('breeding_site','height',','.join([str(h)]*n))
            configuration.config_parser.set('breeding_site','manually_filled',','.join([str(mf)]+[str(0)]*(n-1)))
            model=Model(configuration)
            time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
            utils.showPlot(utils.plot(model,subplots=[{'lwO':'','O':list([90]),'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)), title='Manually Filled:%scm. Height: %scm.(Oct-Nov-Dic just prom available)'%(mf,h))
            #utils.showPlot(utils.plot(model,subplots=[{'E':''}],plot_start_date=datetime.date(2017,10,1)),title='Manually Filled:%scm. Height: %scm.(Oct-Nov-Dic just prom available)'%(mf,h))
            print('mf:%s h:%s Max E: %s'%(mf,h,np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1))))

            #is OEquations perturbing the result somehow?No, the results match.
            #model2=Model(configuration)
            #time_range2,initial_condition2,Y2=model2.solveEquations(method='rk')
            #print(np.linalg.norm((Y[:,:model.parameters.OVIPOSITION.start]-Y2)))

    if(case==2):
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
            mf,h=configuration.getArray('breeding_site','manually_filled')[0],configuration.getArray('breeding_site','height')[0]
            model=Model(configuration)
            time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
            utils.showPlot(utils.plot(model,subplots=[{'E':''}],plot_start_date=datetime.date(2017,10,1)),title='Manually Filled:%scm. Height: %scm.(Oct-Nov-Dic just prom available)'%(mf,h))

try:
    from otero_precipitation_wrapper import ModelWrapper as _Model
except ImportError:
    pass

def runCpp():
    model=_Model('resources/otero_precipitation.cfg')
    Y1=np.array(model.solveEquations())
    print(np.linalg.norm(Y1),Y1.shape)

    model=Model(Configuration('resources/otero_precipitation.cfg'))
    time_range,initial_condition,Y2=model.solveEquations(method='rk' )
    print(np.linalg.norm(Y2),Y2.shape)

    print('||Y1-Y2||=%s'%np.linalg.norm(Y1-Y2))



if(__name__ == '__main__'):
    if(len(sys.argv)>2 and sys.argv[1]=='show'):
            runShow(sys.argv[2])
    elif(len(sys.argv)>2 and sys.argv[1]=='ovishow'):
        runOviShow(sys.argv[2])
    elif(len(sys.argv)>1 and sys.argv[1]=='project'):
        runProject()
    elif(len(sys.argv)>1 and sys.argv[1]=='spatial'):
        runSpatial()
    elif(len(sys.argv)>1 and sys.argv[1]=='cpp'):
        runCpp()
    else:#the default is just a number indicating which test case to run, or none (test case 1 will will be default)
        if(len(sys.argv)<2):
            case=1
        else:
            case=int(sys.argv[1])
        runCases(case)
