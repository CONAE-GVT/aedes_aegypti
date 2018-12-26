import os
import json
import tempfile
import datetime
import numpy as np
import sys
os.chdir('../../../')
sys.path.append('./src')
import utils
from config import Configuration
from otero_precipitation import Model




DOWNLOAD_FOLDER='download/'
SCENARIO_FULLNAMES={'rc':'condiciones_regulares','it':'temperatura_aumentada','dt':'temperatura_disminuida','ip':'precipitacion_aumentada','dp':'precipitacion_disminuida','oc':'contenedores_externos','ic':'contenedores_internos'}
#taken from https://stackoverflow.com/questions/12522844/change-character-set-for-tempfile-namedtemporaryfile
class MyRandomSequence(tempfile._RandomNameSequence):
    characters = "abcdefghijklmnopqrstuvwxyz0123456789"

tempfile._name_sequence = MyRandomSequence()

def createOrderId(prefix='_'):
	file=tempfile.NamedTemporaryFile(dir=DOWNLOAD_FOLDER,prefix=prefix,delete=False)
	file.close()
	return os.path.basename(file.name)

def generateCSV(time_range,columns,names,start_date,end_date,location,scenario):
    augmented_Y=np.zeros((len(time_range),6))#Y will be augmented with weather variables.
    for i,M in enumerate(columns):
        M=np.array(M)
        augmented_Y[:,i]=M[:,1]
    daily_Y=utils.getDailyResults(time_range,augmented_Y,start_date,end_date)

    id=createOrderId()
    filename=location+'_'+start_date.strftime('%Y-%m-%d')+'_'+end_date.strftime('%Y-%m-%d')+ '_'+SCENARIO_FULLNAMES[scenario]+id+'.csv'
    file=open(DOWNLOAD_FOLDER+filename,'w')

    file.write('Fecha,' + ','.join(names)+ '\n')#add the header
    for d,daily_Y_d in enumerate(daily_Y):
        date_d=start_date+datetime.timedelta(days=d)
        file.write(date_d.strftime('%Y-%m-%d')+','+','.join([str(value) for value in daily_Y_d ])+ '\n')
    return filename

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
    #return json.dumps(os.getcwd())
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
                        'weather':[{'name':'Temperatura','data':T},{'name':'Humedad Relativa','data':RH},{'name':'Precipitacion','data':p}],
                        'filename':generateCSV(time_range,[E,L,A,T,RH,p],['Huevos','Larvas','Adultos','Temperatura','Humedad Relativa','Precipitacion'],start_date,end_date,location,scenario)#TODO:find a better way to pass the header
                         })


def download(GET):
    filename=GET.get('filename')
    response=FileResponse(open(DOWNLOAD_FOLDER+filename, 'rb'))
    response['Content-Disposition'] = 'attachment; filename='+'_'.join(filename.split('_')[:-1])+ os.path.splitext(filename)[1]#remove order id
    return response
