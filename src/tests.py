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
from equationsv2 import diff_eqs as diff_eqs_v2
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
    utils.showPlot()

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
    utils.showPlot()

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
    utils.showPlot()

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
def calculateMetrics(V,ovitrap_eggs_i):
    V/=np.max(V)
    ovitrap_eggs_i/=np.nanmax(ovitrap_eggs_i)
    d=np.nansum((ovitrap_eggs_i-V)**2 )/ np.where(~np.isnan(ovitrap_eggs_i),1,0).sum()
    rho,p_value=stats.spearmanr(ovitrap_eggs_i[~np.isnan(ovitrap_eggs_i)], V[~np.isnan(ovitrap_eggs_i)])
    return d,rho,p_value

def runCases(case):
    if(case==1):
        testModel(Configuration('resources/otero_precipitation.cfg'),subplots=[['E'],['A1+A2'], ['W']])
    #Model against ovitraps
    start_date,end_date=utils.getStartEndDates('data/private/ovitrampas_2017-2018.csv')
    configuration=Configuration('resources/otero_precipitation.cfg',
        {'simulation':{
            'start_date':start_date,
            'end_date':end_date
            }
        })
    if(case==2):
        #this test seems to indicate that the lwe produced by EACH ovitrap is independent of the BS_a. Meaning that there's no value in fitting this by BS_a
        BS_d=configuration.getArray('breeding_site','distribution')
        for BS_a in range(20,220,20):
            configuration.config_parser.set('breeding_site','amount',str(BS_a))
            per_ovitrap=lambda lwE: lwE[:,0]/(BS_a*BS_d[0])
            testModel(configuration,subplots=[{'lwE':'','f':[utils.replaceNegativesWithZeros,per_ovitrap]}],title=BS_a)
    if(case==3):
        for alpha0 in np.linspace(1.0,2,10):
            configuration.config_parser.set('biology','alpha0',str(alpha0))
            BS_a=configuration.getFloat('breeding_site','amount')
            BS_d=configuration.getArray('breeding_site','distribution')
            per_ovitrap=lambda lwE: lwE[:,0]/(BS_a*BS_d[0])
            testModel(configuration,subplots=[{'lwE':'','f':[utils.replaceNegativesWithZeros,per_ovitrap]}],title=str(alpha0),figure=False)
    if(case==4):
        configuration=Configuration('resources/otero_precipitation.cfg')#reset config
        for alpha0 in np.linspace(1.0,2,10):
            configuration.config_parser.set('biology','alpha0',str(alpha0))
            model=Model(configuration)
            time_range,Y0,Y=model.solveEquations(equations=utils.MetricsEquations(model,diff_eqs) )
            model=Model(configuration)
            time_range,Y0,Y_rk=model.solveEquations(equations=utils.MetricsEquations(model,diff_eqs),method='rk' )
            model=Model(configuration)
            time_range,Y0,Y_rkf=model.solveEquations(equations=utils.MetricsEquations(model,diff_eqs),method='rkf' )
            print('>'*200)
            print('rk:  e=||Y-Y_rk||=   %s'%(np.linalg.norm(Y-Y_rk)))
            print('rkf: e=||Y-Y_rkf||=  %s'%(np.linalg.norm(Y-Y_rkf)))
            print('     e_rkf/e_rk=     %s'%(np.linalg.norm(Y-Y_rkf)/np.linalg.norm(Y-Y_rk) ) )
            print('<'*200)
    if(case==5):
        per_ovitrap=lambda lwE: lwE[:,0]/(BS_a*BS_d[0]) if(lwE.ndim==2) else lwE
        BS_a=configuration.getFloat('breeding_site','amount')
        BS_d=configuration.getArray('breeding_site','distribution')

        temperatures=np.arange(278,303.5,0.5)
        for j in np.arange(0,1.1,0.1):
            equations.vR_D_298K=np.array([0.24,0.2088,0.384,0.216,0.372])#+j
            equations.vDeltaH_A=np.array([10798.0,26018.0,14931.0,15725.0,15725.0])#+1e5*j
            equations.vDeltaH_H=np.array([100000.0,55990.0,-472379.00,1756481.0,1756481.0]) #-472379 vs. -473379
            equations.vT_1_2H=np.array([14184.0,304.6,148.0,447.2,447.2])
            #for i,label in enumerate(['elr','lpr','par','ovr1(cycle1)','ovr2(cycle2)$\deltaH_A$ %s'%equations.vDeltaH_A[4]]):
                #pl.plot(temperatures,[vR_D(T)[i] for T in temperatures],label=label)

            #pl.plot(temperatures,[vR_D(T)[4] for T in temperatures],label='ovr2(cycle2)')
            pl.plot(temperatures-273.15,[1/vR_D(T)[4] for T in temperatures],label='$ovr2^{-1}$')
            #testModel(configuration,subplots=[['E','A1+A2',[utils.safeAdd]],{'lwE':'','O':[153,154,155],'f':[utils.replaceNegativesWithZeros,per_ovitrap]}],title='',figure=False)
            break
        pl.xlabel('Temperature (K)')
        pl.ylabel('Development Rate ($days^{-1}$)')
        pl.legend(loc=0)

    if(case==6):
        configuration=Configuration('resources/otero_precipitation.cfg')
        model=Model(configuration)
        time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs,method='rk' )
        common(model,'v1')

    if(case==7):
        configuration=Configuration('resources/otero_precipitation.cfg')
        n=len(configuration.getArray('breeding_site','height'))
        configuration.config_parser.set('breeding_site','manually_filled',','.join([str(0.5)] + [str(0)]*(n-1)))
        model=Model(configuration)
        time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs,method='rk' )
        common(model,'v1 ucar mf')

        configuration=Configuration('resources/otero_precipitation.cfg')
        model=Model(configuration)
        time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs,method='rk' )
        common(model,'v1 ucar')

        configuration=Configuration('resources/otero_precipitation.cfg')
        configuration.config_parser.set('location','name','wunderground')
        model=Model(configuration)
        time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs,method='rk' )
        common(model,'v1 wunderground')

        configuration=Configuration('resources/otero_precipitation.cfg')
        configuration.config_parser.set('location','name','wunderground')
        model=getV2Model(Model(configuration))
        time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs_v2,method='rk' )
        common(model,'v2 wunderground')

        for h in [1,5,15,30]:
            configuration=Configuration('resources/otero_precipitation.cfg')
            configuration.config_parser.set('location','name','wunderground')
            n=len(configuration.getArray('breeding_site','height'))
            configuration.config_parser.set('breeding_site','height',','.join([str(h)]*n))
            configuration.config_parser.set('location','name','wunderground')
            model=getV2Model(Model(configuration))
            time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs_v2,method='rk' )
            common(model,'v2 wunderground h=%s'%h)

    if(case==8):
        configuration=Configuration('resources/otero_precipitation.cfg')#Configuration('resources/1c.cfg')
        configuration.config_parser.set('location','name','wunderground')
        model=getV2Model(Model(configuration))
        time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs_v2,method='rk' )

        error_E0,error_lwE0,error_OV0=[],[],[]
        for ovitrap_id in range(1,151):
            OVITRAP_FILENAME='data/private/ovitrampas_2017-2018.csv'
            ovitrap_eggs_i=utils.getOvitrapEggsFromCsv(OVITRAP_FILENAME,start_date,end_date,ovitrap_id)
            ovitrap_eggs_i=equation_fitter.populate(model.time_range,ovitrap_eggs_i)
            ovitrap_eggs_i=np.array(ovitrap_eggs_i,dtype=np.float)#this change None for np.nan

            BS_l=model.parameters.BS_l
            E0=Y[:,0:BS_l].sum(axis=1)
            lwE0=np.array([Y[(np.abs(time_range-t)).argmin(),BS_l]-Y[(np.abs(time_range-(t-7))).argmin(),BS_l] for t in time_range])

            OV0=[]
            vBS_d=model.parameters.vBS_d
            for i,t in enumerate(time_range):
                T_t=model.parameters.weather.T(t)
                vW_t=model.parameters.vW(t)
                egn,ovr1,ovr2,A1,A2=63.,equations.vR_D(T_t)[0],equations.vR_D(T_t)[0],Y[i,-2],Y[i,-1]
                OV_t= egn*( ovr1 *A1  + ovr2* A2)*equations.f(vW_t,vBS_d)
                OV0=OV0 + [OV_t[0]]
            OV0=np.array(OV0)

            square,rho,p_value=calculateMetrics(E0,ovitrap_eggs_i)
            error_E0=error_E0+[[square,rho,p_value]]
            square,rho,p_value=calculateMetrics(lwE0,ovitrap_eggs_i)
            error_lwE0=error_lwE0+[[square,rho,p_value]]
            square,rho,p_value=calculateMetrics(OV0,ovitrap_eggs_i)
            error_OV0=error_OV0+[[square,rho,p_value]]
            #print('ovitrap %s Error: %s rho: %s p-value: %s'%(ovitrap_id,error,rho,p_value) )

        for error in [error_E0,error_lwE0,error_OV0]:
            pl.figure()
            error=np.array(error)
            error[:,0]=error[:,0]/error[:,0].max()
            for i,label in enumerate(['square','rho','p_value']):
                pl.plot(range(1,151),error[:,i],label=label)
                pl.legend(loc=0)

    if (case==9):
        classes=[line.strip().split(',') for line in open('data/private/dtw_results.csv','r').readlines()]
        for method in [1,2]:
            pl.figure()
            for class_number in [1,2,3]:
                pl.title('method %s class%s, (1:dtw,2:kmeans)'%(method,class_number))
                pl.subplot(300 + 10 + class_number)
                for ovitrap_id in range(1,151):
                    OVITRAP_FILENAME='data/private/ovitrampas_2017-2018.csv'
                    start_date,end_date=utils.getStartEndDates(OVITRAP_FILENAME)
                    ovitrap_eggs_i=utils.getOvitrapEggsFromCsv(OVITRAP_FILENAME,start_date,end_date,ovitrap_id)
                    #take away Nones
                    ovitrap_eggs=[e for e in ovitrap_eggs_i if e!=None]
                    ovitrap_days=[datetime.timedelta(days=d)+datetime.datetime.combine(start_date,datetime.time()) for d in range(0,len(ovitrap_eggs_i)) if ovitrap_eggs_i[d]!=None]
                    if(int(classes[ovitrap_id][method])==class_number):
                        pl.plot(ovitrap_days, ovitrap_eggs, '-')

    if (case==10):
        vStd_accum=[1e10]*152
        for ovitrap_id in range(1,151):
            OVITRAP_FILENAME='data/private/ovitrampas_2017-2018.full.csv'
            values=utils.getOvitrapEggsFromCsv2(OVITRAP_FILENAME,None,None,ovitrap_id)
            ovitrap_days=values.keys()
            ovi_a=[values[date][0] for date in ovitrap_days]
            ovi_b=[values[date][1] if len(values[date])>1 else values[date][0] for date in ovitrap_days]
            ovi_std=[np.std([ovi_a[i],ovi_b[i]]) if (ovi_a[i]!= None and ovi_b[i]!= None) else None for i in range(0,len(ovi_a))]
            #pl.plot(ovitrap_days, ovi_a, '-')
            #pl.plot(ovitrap_days, ovi_b, '-')
            pl.plot(ovitrap_days, ovi_std, '-')
            ovi_std_safe=[ovi_std[i] for i in range(len(ovitrap_days)) if ovi_a[i] is not None and ovi_b[i] is not None]
            if(len(ovi_std_safe) != len(ovi_std)): continue
            vStd_accum[ovitrap_id]=sum(ovi_std_safe);


        sorted_ovi=np.argsort(vStd_accum)
        for ovitrap_id in sorted_ovi:
            print('Ovitrap id:%s, Accumulated std:%s'%(ovitrap_id,vStd_accum[ovitrap_id]))
            #TODO:find the one with min std and max std, and compare those to model.


    utils.showPlot()

def getV2Model(model):
    #modify some parameters to make them compatible with v2
    configuration=model.configuration
    parameters=model.parameters
    n=parameters.n
    parameters.BS_l=int(configuration.getFloat('breeding_site','levels'))#BS levels#TODO:by changing this value, we get a different number of adults, which looks too big. Check if there isn't an error somewhere, or a way to make it more stable
    BS_l=parameters.BS_l
    parameters.EGG=slice(0,n*BS_l)#in R^n
    parameters.LARVAE=slice(n*BS_l,(1+BS_l)*n)#in R^n
    parameters.PUPAE=slice((1+BS_l)*n,(2+BS_l)*n)#in R^n
    parameters.ADULT1=(2+BS_l)*n#in R
    parameters.ADULT2=(2+BS_l)*n+1#in R
    #parameters.TODO:add the B matrix helper
    parameters.mBS_l=np.repeat(range(0,BS_l),n).reshape((BS_l,n))
    E0v2=parameters.initial_condition[:n].repeat(BS_l)/BS_l
    parameters.initial_condition=np.insert(parameters.initial_condition[n:],0,E0v2)

    #experimental
    p=configuration.getFloat('simulation','egn_corrector_p')#TODO:change p for something with meaning...
    parameters.egnCorrector=utils.EgnCorrector(p,model.parameters.BS_a,model.start_date,model.end_date)
    return model

def runModelv2(case):
    configuration=Configuration('resources/otero_precipitation.cfg')
    model=getV2Model(Model(configuration))

    if(case==1):
        for dis in [0.9,0.5,0.1]:#this test actually shows something meaningfull when you have just 2 types of identical BS, one with mf and one without.
            n=model.parameters.n
            model.parameters.vBS_d=np.array([dis]+ [(1-dis)/(n-1)]*(n-1))
            time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs_v2,method='rk' )
            utils.plot(model,subplots=[{'E':'','O':[153],'f':[utils.safeAdd,utils.replaceNegativesWithZeros,utils.safeNormalize]},['A1+A2',[utils.safeNormalize] ],['W'],['p']],title='BS_d: %s'%model.parameters.vBS_d.round(2))
            print('egn Correction mean: %s'%np.mean(parameters.egnCorrector.egn_corrections))

    if(case==6):
        time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs_v2,method='rk' )
        common(model,'v2')

    if(case==7):
        configuration=Configuration('resources/1c.cfg')#use this config.
        model=getV2Model(Model(configuration))
        time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs_v2,method='rk' )
        utils.plot(model,subplots=[['E'],['W']])

    for warning in model.warnings:
        print('# WARNING: ' + warning)
    utils.showPlot()

def common(model,title):
    print('Max E: %s'%np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1)))
    utils.plot(model,subplots=[
        {'E':'',    'O':[153],'f':[utils.safeAdd,utils.replaceNegativesWithZeros,utils.safeNormalize]},
        {'A1+A2':'','O':[153],'f':[utils.safeAdd,utils.replaceNegativesWithZeros,utils.safeNormalize]},
        ['lwE',[utils.safeAdd,utils.replaceNegativesWithZeros,utils.safeNormalize]]  ,['W']],title=title)#, ['W'],['p'],['T']

if(__name__ == '__main__'):
    if(len(sys.argv)>2 and sys.argv[1]=='show'):
            runShow(sys.argv[2])
    elif(len(sys.argv)>2 and sys.argv[1]=='ovishow'):
        runOviShow(sys.argv[2])
    elif(len(sys.argv)>1 and sys.argv[1]=='project'):
        runProject()
    elif(len(sys.argv)>1 and sys.argv[1]=='spatial'):
        runSpatial()
    elif(len(sys.argv)>1 and sys.argv[1]=='v2'):
        runModelv2(int(sys.argv[2]))
    else:#the default is just a number indicating which test case to run, or none (test case 1 will will be default)
        if(len(sys.argv)<2):
            case=1
        else:
            case=int(sys.argv[1])
        runCases(case)
