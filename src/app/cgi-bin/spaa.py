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
    configuration=Configuration('resources/otero_precipitation.cfg',
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

    time_range,initial_condition,Y= model.solveEquations(method='rk')
    EGG,LARVAE,ADULT1,ADULT2=model.parameters.EGG,model.parameters.LARVAE,model.parameters.ADULT1,model.parameters.ADULT2
    BS_a=configuration.getFloat('breeding_site','amount')

    E=[ [ daysSinceEpoch(start_datetime,t), np.sum(Y[i,EGG])/BS_a ] for i,t in enumerate(time_range)]
    L=[ [ daysSinceEpoch(start_datetime,t), np.sum(Y[i,LARVAE])/BS_a ] for i,t in enumerate(time_range)]
    A=[ [ daysSinceEpoch(start_datetime,t), (Y[i,ADULT1]+Y[i,ADULT2])/BS_a ] for i,t in enumerate(time_range)]

    weather=model.parameters.weather
    p=[ [ daysSinceEpoch(start_datetime,t), weather.p(t) ] for t in time_range]
    T=[ [ daysSinceEpoch(start_datetime,t), np.asscalar(weather.T(t)) - 273.15 ] for t in time_range]
    RH=[ [ daysSinceEpoch(start_datetime,t), np.asscalar(weather.RH(t)) ] for t in time_range]

    return json.dumps({
                        'population':[{'name':'Huevos','data':E},{'name':'Larvas','data':L},{'name':'Adultos','data':A}],
                        'weather':[{'name':'Temperatura','data':T},{'name':'Humedad Relativa','data':RH},{'name':'Precipitacion','data':p}]
                        })
