import os
import json
import datetime
import numpy as np
import sys
os.chdir('../../../')
sys.path.append('./src')
import utils
from config import Configuration
from otero_precipitation import Model
from equations import diff_eqs

SCENARIO_FULLNAMES={'rc':'condiciones_regulares','it':'temperatura_aumentada','dt':'temperatura_disminuida','ip':'precipitacion_aumentada','dp':'precipitacion_disminuida','oc':'contenedores_externos','ic':'contenedores_internos'}

def daysSinceEpoch(start_datetime,days):
    epoch = datetime.datetime.utcfromtimestamp(0)
    return (start_datetime + datetime.timedelta(days=days) - epoch).total_seconds() * 1000.0

def applyScenario(model,scenario):
    if(scenario in ['it','dt']):
        if(scenario=='it'):#increased temperature
            a=1.2
        elif(scenario=='dt'):#decreased temperature
            a=0.8
        T=model.parameters.weather.T
        model.parameters.weather.T=lambda t: (T(t)-273.15)*a + 273.15

    if(scenario in ['ip','dp']):
        if(scenario=='ip'):#increased precipitation
            a=1.2
        elif(scenario=='dp'):#decreased precipitation
            a=0.8
        p=model.parameters.weather.p
        model.parameters.weather.p=lambda t: p(t)*a


    return model

def getConfig(configuration,scenario):
    if(scenario in ['oc','ic']):
        if(scenario=='oc'):#80% outside containers
            a=0.8
        elif(scenario=='ic'):#80% inside containers
            a=0.2
        vBS_od=configuration.getArray('breeding_site','outside_distribution')
        vBS_id=configuration.getArray('breeding_site','inside_distribution')

        vBS_od=a*vBS_od/np.sum(vBS_od)
        vBS_id=(1-a)*vBS_id/np.sum(vBS_id)
        assert 1-(np.sum(vBS_od)+np.sum(vBS_id))<1e-6, 'vBS_od+vBS_id=%s'%(np.sum(vBS_od)+np.sum(vBS_id))
        configuration.config_parser.set('breeding_site','outside_distribution',','.join([str(value) for value in vBS_od ]))
        configuration.config_parser.set('breeding_site','inside_distribution',','.join([str(value) for value in vBS_id ]))

    return configuration

def runSimulation(GET):
    start_datetime=datetime.datetime.strptime(GET.get('start_date'),'%Y-%m-%d')
    start_date=start_datetime.date()
    end_date=datetime.datetime.strptime(GET.get('end_date'),'%Y-%m-%d').date()
    location=GET.get('location')
    scenario=GET.get('scenario')
    configuration=Configuration('resources/1c.cfg',
        {'simulation':{
            'start_date':start_date,
            'end_date':end_date,
            },
        'location':{
            'name':str(location)+'.full'#unicode to str
        },
        })
    configuration=getConfig(configuration,scenario)
    model=applyScenario(Model(configuration), scenario)

    time_range,initial_condition,Y= model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
    EGG,LARVAE,ADULT1,ADULT2,OVIPOSITION=model.parameters.EGG,model.parameters.LARVAE,model.parameters.ADULT1,model.parameters.ADULT2,model.parameters.OVIPOSITION
    BS_a=configuration.getFloat('breeding_site','amount')
    indexOf=lambda t: (np.abs(time_range-t)).argmin()

    E=[ [ daysSinceEpoch(start_datetime,t), np.sum(Y[i,EGG])/BS_a ] for i,t in enumerate(time_range)]
    L=[ [ daysSinceEpoch(start_datetime,t), np.sum(Y[i,LARVAE])/BS_a ] for i,t in enumerate(time_range)]
    A=[ [ daysSinceEpoch(start_datetime,t), (Y[i,ADULT1]+Y[i,ADULT2])/BS_a ] for i,t in enumerate(time_range)]

    lwO=np.array([Y[indexOf(t),OVIPOSITION]-Y[indexOf(t-7),OVIPOSITION] for t in time_range])
    lwO_mean=np.array([lwO[indexOf(t-7):indexOf(t+7)].mean(axis=0) for t in time_range])
    O=[ [ daysSinceEpoch(start_datetime,t), np.sum(lwO_mean[i])/BS_a ] for i,t in enumerate(time_range)]

    weather=model.parameters.weather
    precipitations = utils.getPrecipitationsFromCsv('data/public/'+location+'.full.csv',start_date,end_date)
    p=[ [ daysSinceEpoch(start_datetime,t), precipitations[int(t)] ] for t in time_range]
    T=[ [ daysSinceEpoch(start_datetime,t), np.asscalar(weather.T(t)) - 273.15 ] for t in time_range]
    RH=[ [ daysSinceEpoch(start_datetime,t), np.asscalar(weather.RH(t)) ] for t in time_range]

    return json.dumps({
                        'population':[{'name':'Huevos','data':E,'type':'scatter'},{'name':'Larvas','data':L,'type':'scatter'},{'name':'Adultos','data':A,'type':'scatter'},{'name':'Oviposicion','data':O,'type':'scatter'}],
                        'weather':[{'name':'Temperatura','data':T,'type':'scatter'},{'name':'Humedad Relativa','data':RH,'type':'scatter'},{'name':'Precipitacion','data':p,'type':'bar'}]
                        })
