#coding: utf-8
import os
import re
import sys
import utils
import datetime
import numpy as np
from config import Configuration
from configparser import ConfigParser
from otero_precipitation import Model
from equations import diff_eqs,vR_D
import equations
from spatial_equations import diff_eqs as spatial_diff_eqs
import pylab as pl
import similaritymeasures as sm
import plotly.graph_objs as go
from plotly import tools

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

from scipy import stats
def calculateMetrics(time_range,lwO_mean,ovitrap_eggs_i):
    valid_ovi_idx=~np.isnan(ovitrap_eggs_i)
    rmse=utils.rmse(ovitrap_eggs_i[valid_ovi_idx],lwO_mean[valid_ovi_idx])
    cort=utils.cort(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx])
    pearson,p_value=stats.pearsonr(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx])

    reversed_valid_ovi_idx=valid_ovi_idx[::-1]
    first,last=np.argmax(valid_ovi_idx), len(reversed_valid_ovi_idx)-np.argmax(reversed_valid_ovi_idx)-1
    x=np.array([[time_range[idx],lwO_mean[idx]] for idx in range(first,last)])
    y=np.array([ [time_range[idx],ovitrap_eggs_i[idx] ] for idx,isValid in enumerate(valid_ovi_idx) if isValid])
    fd=sm.frechet_dist(x,y)
    dtw, path = sm.dtw(x, y)
    D_1=utils.D(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx],k=1)
    D=utils.D(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx])
    D_2=utils.D(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx],k=2)
    D_4=utils.D(ovitrap_eggs_i[valid_ovi_idx], lwO_mean[valid_ovi_idx],k=4)

    return rmse, cort,pearson,fd,dtw,D,D_1,D_2,D_4

import equation_fitter
def runCases(case):
    if(case==0):
        ovi_range=range(1,151)
        errors_by_height=np.empty((15,151,9))
        errors_by_height[:]=np.nan
        for h in range(1,15):
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')#TODO:fix data and
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(20)))# uncomment these two
            configuration.config_parser.set('breeding_site','height',str(h))
            model=Model(configuration)
            time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')

            #errors=[[1e15,-1,-1,1e15,1e15,1e15]]*151#just to fill the ovitrap 0 that do not exist in reality
            for ovitrap_id in ovi_range:
                OVITRAP_FILENAME='data/private/ovitrampas_2017-2018.full.csv'
                values=utils.getOvitrapEggsFromCsv2(OVITRAP_FILENAME,None,None,ovitrap_id)
                ovitrap_days=values.keys()
                dates=[model.start_date + datetime.timedelta(t) for t in time_range]

                ovi=[utils.noneMean(values[date]) if date in values else None for date in dates]#TODO:WARNING!this will repeat values if model granularity is not 1 value per day.
                ovi=np.array(equation_fitter.populate(model.time_range,ovi))
                ovi=np.array(ovi,dtype=np.float)#this change None for np.nan

                indexOf=lambda t: (np.abs(time_range-t)).argmin()
                OVIPOSITION=model.parameters.OVIPOSITION
                BS_a=model.parameters.BS_a
                O=Y[:,OVIPOSITION]
                lwO=np.sum([Y[indexOf(t),OVIPOSITION]-Y[indexOf(t-7),OVIPOSITION] for t in time_range],axis=1)/BS_a#calculate the difference,sum up all levels, and divide by amount of containers

                errors_by_height[h][ovitrap_id]=calculateMetrics(time_range,lwO,ovi)


            errors=errors_by_height[h]
            rmse, cort,pearson,fd,dtw,D,D_1,D_2,D_4=[errors[:,i] for i in range(errors.shape[1])]
            print('''h: %scm.
                         id,    score
                rmse:    %3s,   %s
                cort:    %3s,   %s
                pearson: %3s,   %s
                fd:      %3s,   %s
                dtw:     %3s,   %s
                D:       %3s,   %s'''%
                (h,
                np.nanargmin(rmse),np.nanmin(rmse),
                np.nanargmin(cort),np.nanmin(cort),
                np.nanargmin(pearson),np.nanmin(pearson),
                np.nanargmin(fd),np.nanmin(fd),
                np.nanargmin(dtw),np.nanmin(dtw),
                np.nanargmin(D),np.nanmin(D)
                ) )

            #print first N
            for i in range(errors.shape[1]):
                f=errors[:,i]
                ovi_sorted=np.argsort(f)
                N=5
                print('''h: %scm.
                     id,    score
                     f:    %3s  <---->  %s'''%
                     (h,
                     ', '.join(map(str,ovi_sorted[:N])), ', '.join(map(str,f[ovi_sorted[:N]]))
                     ))
                #print last N(excluding the ficticius one)
                print('''h: %scm.
                     id,    score
                     f:    %3s  <---->  %s'''%
                     (h,
                     ', '.join(map(str,ovi_sorted[-N-1:-1])), ', '.join(map(str,f[ovi_sorted[-N-1:-1]]))
                     ))

        np.save('errors_by_height.npy',errors_by_height)
        print(errors_by_height.shape)
        pl.show()

    if(case==1):
        h=1.
        configuration=Configuration('resources/1c.cfg')
        configuration.config_parser.set('location','name','cordoba.full')
        configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
        configuration.config_parser.set('breeding_site','height',str(h))
        model=Model(configuration)
        time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
        utils.showPlot(utils.plot(model,subplots=[{'cd':'','lwO':'','O':list([34,19,133,1,56,16,25,143,59,44]),'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
        title='Height: %scm.'%h,
        xaxis_title='Fecha',
        yaxis_title='Nº de huevos')

        #utils.showPlot(utils.plot(model,subplots=[{'E':''}],plot_start_date=datetime.date(2017,10,1)),title='Manually Filled:%scm. Height: %scm.(Oct-Nov-Dic just prom available)'%(mf,h))
        #utils.showPlot(utils.plot(model,subplots=[{'pa':''}]))
        print('h:%s Max E: %s'%(h,np.max(np.sum(model.Y[:,model.parameters.EGG],axis=1))))
        print(model.warnings)

        #is OEquations perturbing the result somehow?No, the results match.
        #model2=Model(configuration)
        #time_range2,initial_condition2,Y2=model2.solveEquations(method='rk')
        #print(np.linalg.norm((Y[:,:model.parameters.OVIPOSITION.start]-Y2)))

    if(case==2):
        for mf  in [0.,3.]:
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
            configuration.config_parser.set('breeding_site','manually_filled',str(mf))
            h=configuration.getArray('breeding_site','height')[0]
            model=Model(configuration)
            time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
            utils.showPlot(utils.plot(model,subplots=[{'E':''}],plot_start_date=datetime.date(2017,10,1)),title='Manually Filled:%scm. Height: %scm.'%(mf,h))
            print(model.warnings)

    if(case==3):
        for year in range(2016,2019):
            for month in range(1,13):
                start_date,end_date=datetime.date(year,month,1),datetime.date(year+int(month/12),month%12 +1,1)
                precipitations = utils.getPrecipitationsFromCsv(sys.argv[2],start_date,end_date)
                print('period: %s to %s        %smm.'%(start_date,end_date,np.sum(precipitations)))
    if(case==4):
        ovi_range=range(1,151)
        ovi_mean=[1e10]*151
        for ovitrap_id in ovi_range:
            OVITRAP_FILENAME='data/private/ovitrampas_2017-2018.full.csv'
            values=utils.getOvitrapEggsFromCsv2(OVITRAP_FILENAME,None,None,ovitrap_id)
            dates=values.keys()
            ovi_a=[values[date][0] if date in values else None for date in dates]#TODO:WARNING!this will repeat values if model granularity is not 1 value per day.
            ovi_a=np.array(ovi_a,dtype=np.float)#this change None for np.nan
            ovi_b=[values[date][1] if (date in values and len(values[date])>1) else None for date in dates]
            ovi_b=np.array(ovi_b,dtype=np.float)#this change None for np.nan
            ovi_mean[ovitrap_id]=np.nanmean(np.abs(ovi_a-ovi_b)/(ovi_a+ovi_b))

        ovi_mean=np.array(ovi_mean)
        ovi_ordered=np.argsort(ovi_mean)
        print(ovi_ordered)
        print(ovi_mean[ovi_ordered])
    if(case==5):
        for location in ['cordoba.full.weather-2019-01-01','cordoba.full.weather-2019-02-01']:
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name',location)
            start_date,end_date=utils.getStartEndDates('data/public/'+location+'.csv')
            configuration.config_parser.set('simulation','end_date',str(end_date))
            model=Model(configuration)
            time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
            utils.showPlot(utils.plot(model,subplots=[{'cd':'','lwO':'','O':list([143]),'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
            title=location.replace('.full.weather-',' ').title(),
            xaxis_title='Date',
            yaxis_title='Number of eggs')
            print(model.warnings)

    if(case==6):
        config_parser = ConfigParser()
        config_parser.read('resources/get_weather.cfg')
        for location in config_parser.sections():
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name',location+'.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
            model=Model(configuration)
            time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
            utils.showPlot(utils.plot(model,subplots=[{'lwO':'','f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
            title='%s '%(location.replace('_',' ').title()),
            xaxis_title='Date',
            yaxis_title='Number of eggs')
            print(model.warnings)

    if(case==7):
        configuration=Configuration('resources/1c.cfg')
        configuration.config_parser.set('location','name','cordoba.full')
        configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
        model=Model(configuration)
        time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
        utils.showPlot(utils.plot(model,subplots=[{'E':':','A1+A2':'','f':[utils.safeAdd]}]),
            title='Cordoba',
            xaxis_title='Date',
            yaxis_title='Individuals')
        print(model.warnings)

    if(case==8):
        LOCATION='cordoba'
        PLOT_START_DATE=datetime.date(2018,12,15)#maybe by parameter?
        FORECAST=51
        PLOT_END_DATE=PLOT_START_DATE+datetime.timedelta(FORECAST)
        DATA_FOLDER='data/public/'
        HISTORY_FOLDER=DATA_FOLDER  + '.history/'
        data_A,data_O,data_W=[],[],[]
        filenames=os.listdir(HISTORY_FOLDER)
        filenames.sort()
        i=-1
        for filename in  filenames:
            if(filename=='.empty' or not filename.startswith(LOCATION)): continue
            location,year,month,day=filename.replace('.full.weather','').replace('.csv','').split('-')
            simulation_date=datetime.date(int(year),int(month),int(day))
            if( not (PLOT_START_DATE<=simulation_date<=PLOT_END_DATE)): continue#the equals is because we want to  have one last curve with no forecast
            i=i+1
            if(i==0):
                color='rgb(0, 0, 0)'
            elif(i==FORECAST-1):
                color='rgb(255, 0, 255)'
            else:
                color = 'rgb(%s, %s, 0)'%(int( (1-i/FORECAST) * 255),int(i/FORECAST * 255))
            if(i%7!=0 and i!=FORECAST-1): continue#every 7 days but the last one must be shown.
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','.history/'+filename.replace('.csv',''))
            configuration.config_parser.set('simulation','end_date',str(PLOT_END_DATE))
            model=Model(configuration)
            time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
            data_A+=utils.plot(model,subplots=[{'A1+A2':str(simulation_date),'f':[utils.safeAdd]}],plot_start_date=PLOT_START_DATE,color=color)
            data_O+=utils.plot(model,subplots=[{'lwO':str(simulation_date)  ,'f':[utils.safeAdd]}],plot_start_date=PLOT_START_DATE,color=color)
            data_W+=utils.plot(model,subplots=[{'W':str(simulation_date)    ,'f':[]             }],plot_start_date=PLOT_START_DATE,color=color)
            print(filename,model.warnings)


        x=[]
        y=[]
        for serie in data_A:
            x+= [( datetime.datetime.strptime(data_A[-1]['name'],'%Y-%m-%d') - datetime.datetime.strptime(serie['name'],'%Y-%m-%d') ).days]#TODO:we depend on using simulation date as name
            S_last,S_i=np.array(data_A[-1]['y']),np.array(serie['y'])
            y+= [np.sum(np.abs(S_last-S_i)/S_last)*1/len(S_i) ]#relative difference mean. 0<=||S_last-S_i||/(S_last+S_i) <= 1 by triangular inequality. => sum(...)/n in [0,1].
            assert len(S_i)==len(S_last)
        x,y=[e for e in reversed(x)],[e for e in reversed(y)]

        utils.showPlot(data_A,title='Adults in '+LOCATION.title(),xaxis_title='Date',yaxis_title='Individuals')
        utils.showPlot(data_O,title='Oviposition in '+LOCATION.title(),xaxis_title='Date',yaxis_title='Eggs')
        utils.showPlot(data_W,title='Water in '+LOCATION.title(),xaxis_title='Date',yaxis_title='cm.')
        utils.showPlot([go.Scatter(x=x,y=y, name='')],title='',xaxis_title='Amount of days forecast',yaxis_title='Mean relative difference(?)')

    names=['rmse', 'cort','pearson','fd','dtw','D','D_1','D_2','D_4']
    if(case==9):
        errors_by_height=np.load('errors_by_height.npy')
        for d,name in enumerate(names):
            utils.showPlot([go.Surface(z=errors_by_height[:,:,d] )],title=name,scene=dict(xaxis=dict(title='Ovitrap id'),yaxis=dict(title='Height')) )
    if(case==10):
        errors_by_height=np.load('errors_by_height.npy')
        d=5
        for h in range(1,15):
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(30)))
            configuration.config_parser.set('breeding_site','height',str(h))
            model=Model(configuration)
            time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
            o_h=[]
            for i in range(1,151):
                if np.nanargmin(errors_by_height[:,i,d])==h: o_h+=[i]

            if(o_h):
                utils.showPlot(utils.plot(model,subplots=[{'cd':'','lwO':'','O':list(o_h),'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1)),
                    title=names[d]+' Height: %scm.'%h,
                    xaxis_title='Fecha',
                    yaxis_title='Nº de huevos')

    if(case==11):
        errors_by_height=np.load('errors_by_height.npy')
        for d,name in enumerate(names):
            fig = tools.make_subplots(rows=errors_by_height.shape[0], cols=1)
            for h in range(1,errors_by_height.shape[0]):
                fig.append_trace(go.Scatter(x=np.array(range(0,errors_by_height.shape[1])), y=errors_by_height[h,:,d],name='%scm.'%h),h,1)
                fig['layout']['yaxis'+str(h)].update(range=[1,np.nanmax(errors_by_height[:,:,d])])
            fig['layout']['title']=name
            utils.showPlot(fig)

    if(case==12):
        errors_by_height=np.load('errors_by_height.npy')
        d=5#this has to be the same as test 10!
        fig = tools.make_subplots(rows=2, cols=2)
        for i,h in enumerate([1,4,6,8]):
            configuration=Configuration('resources/1c.cfg')
            configuration.config_parser.set('location','name','cordoba.full')
            configuration.config_parser.set('simulation','end_date',str(datetime.date.today()+datetime.timedelta(30)))
            configuration.config_parser.set('breeding_site','height',str(h))
            model=Model(configuration)
            time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
            fig.append_trace(utils.plot(model,subplots=[{'cd':'','lwO':'','O':[np.nanargmin(errors_by_height[h,:,d])],'f':[utils.safeAdd]}],plot_start_date=datetime.date(2017,10,1))[0],int(i/2) +1,i%2 +1)
            #fig['layout']['yaxis'+str(h)].update(range=[1,np.nanmax(errors_by_height[:,:,d])])
        fig['layout']['title']=''
        utils.showPlot(fig)


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

import netCDF4 as nc
def runInfo(nc_filename):
    grp = nc.Dataset(nc_filename)
    lats = grp.variables['lat'][:]
    lons = grp.variables['lon'][:]
    precipitations=grp.variables['precipitationCal']
    lat=-31.420082
    lon=-64.188774
    p=precipitations[(abs(lons-lon)).argmin(),(abs(lats-lat)).argmin()]
    print(precipitations.shape)
    print(grp)
    print(p)

def rewritehistory():
    DATA_FOLDER='data/public/'
    HISTORY_FOLDER=DATA_FOLDER  + '.history/'
    for filename in  os.listdir(HISTORY_FOLDER):
        if(filename=='.empty' or not open(HISTORY_FOLDER+filename).readline().startswith('Date')): continue
        location,year,month,day=filename.replace('.full.weather','').replace('.csv','').split('-')
        start_date,tmp=utils.getStartEndDates(HISTORY_FOLDER+filename)
        end_date=datetime.date(int(year),int(month),int(day))
        precipitations=utils.getPrecipitationsFromCsv(DATA_FOLDER+location+'.full.csv',start_date,end_date)
        content='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h\n'
        for i,line in enumerate(open(HISTORY_FOLDER+filename).readlines()[1:]):
            fields=line.rstrip().split(',')
            if(i<len(precipitations)): fields[4]=str(precipitations[i])
            content+=','.join(fields)+',,\n'
        open('data/public/out/'+filename,'w').write(content)
        #print(start_date,end_date)


if(__name__ == '__main__'):
    if(len(sys.argv)>1 and sys.argv[1]=='spatial'):
        runSpatial()
    elif(len(sys.argv)>1 and sys.argv[1]=='cpp'):
        runCpp()
    elif(len(sys.argv)>1 and sys.argv[1]=='info'):
        runInfo(sys.argv[2])
    elif(len(sys.argv)>1 and sys.argv[1]=='rewrite'):
        rewritehistory()
    else:#the default is just a number indicating which test case to run, or none (test case 1 will will be default)
        if(len(sys.argv)<2):
            case=1
        else:
            case=int(sys.argv[1])
        runCases(case)
