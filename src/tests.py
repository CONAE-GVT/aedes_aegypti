#coding: utf-8
import os
import re
import sys
import utils
import datetime
import numpy as np
from config import Configuration
from otero_precipitation import Model
from equations import diff_eqs
from spatial_equations import diff_eqs as spatial_diff_eqs

def testModel(configuration, subplots=[],plot_start_date=None,title='',figure=True,color=None):
    model=Model(configuration)
    model.parameters.calls=None
    model.parameters.negatives=None
    time_range,INPUT,RES=model.solveEquations(equations=utils.MetricsEquations(model,diff_eqs),method='rk' )
    utils.plot(model,subplots,plot_start_date,title,figure,color)

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
    parameters.P=utils.getPreferenceMatrix()
    HEIGHT,WIDTH=parameters.P.shape[:2]
    tmp=np.zeros((HEIGHT,WIDTH,3*n + 3))
    tmp[int(HEIGHT/2),int(WIDTH/2),:]=1.
    parameters.initial_condition=np.append(parameters.initial_condition,[0])#append flyers#TODO:use a config file
    parameters.initial_condition=(parameters.initial_condition*tmp).reshape((HEIGHT*WIDTH*(3*n + 3) ))#TODO:check that this does what we expect.
    parameters.vBS_a=parameters.BS_a*np.ones((HEIGHT,WIDTH))#np.random.random((WIDTH,HEIGHT))#TODO:do something about this...
    parameters.vBS_d=parameters.vBS_d*np.ones((HEIGHT,WIDTH,n))
    parameters.vAlpha=(parameters.vAlpha0*np.ones((HEIGHT,WIDTH,n)) )/parameters.vBS_a[:,:,np.newaxis]
    parameters.ovr=np.ones((HEIGHT,WIDTH))/0.229#TODO:implement!!!!

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

    utils.showPlot()

if(__name__ == '__main__'):
    if(len(sys.argv)>2 and sys.argv[1]=='show'):
            runShow(sys.argv[2])
    elif(len(sys.argv)>2 and sys.argv[1]=='ovishow'):
        runOviShow(sys.argv[2])
    elif(len(sys.argv)>1 and sys.argv[1]=='project'):
        runProject()
    elif(len(sys.argv)>1 and sys.argv[1]=='spatial'):
        runSpatial()
    else:#the default is just a number indicating which test case to run, or none (test case 1 will will be default)
        if(len(sys.argv)<2):
            case=1
        else:
            case=int(sys.argv[1])
        runCases(case)
