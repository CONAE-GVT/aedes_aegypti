import os
import re
import sys
import logging
import datetime
import urllib2
import gdas_lib
import imerg_lib
import pygrib
import numpy as np
import netCDF4 as nc
from utils import daterange,getLocations
from configparser import ConfigParser
import multiprocessing as mp

DATA_FOLDER='data/public/'
IMERG_FOLDER=DATA_FOLDER+'/imerg/'
GDAS_FOLDER=DATA_FOLDER+'/gdas/'
FORECAST_FOLDER=DATA_FOLDER+'/forecast/'
FORECAST_TRH_FOLDER=FORECAST_FOLDER+'/T_RH/'
FORECAST_P_FOLDER=FORECAST_FOLDER+'/P/'
LOG_FILENAME='logs/get_weather.log'
logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s',filename=LOG_FILENAME,level=logging.DEBUG)

def renameAll(dsts):
    for dst in dsts:
        for filename in os.listdir(dst):
            if(not os.path.isfile(dst+'/'+filename) or filename.startswith('gdas1')): continue
            date_str,hours_str=re.findall(r'GFS_([0-9]{8})00\+([0-9]{3}).grib2',filename)[0]
            a_date_time=datetime.datetime.strptime(date_str,'%Y%m%d') + datetime.timedelta(hours=int(hours_str))
            new_filename=gdas_lib.getFilename(a_date_time.date(),'%02d'%a_date_time.hour)
            os.rename(dst+'/'+filename,dst+'/'+new_filename)

def downloadAll(url,folder):
    html_list=urllib2.urlopen(url).read()
    filenames=re.findall(r'.*\<a href=\"([A-Za-z0-9\+\._]+)\"\>.*',html_list)
    for filename in filenames:
        if(filename.endswith('.grib2_NACC_2m') or filename.endswith('.grib2_ACC')):
            open(folder+'/'+filename, 'w').write(urllib2.urlopen(url+'/'+filename).read())

def downloadForecast():
    #download
    today=datetime.date.today()
    src_T_RH = 'http://meteo.caearte.conae.gov.ar/GFS_025/DATA/%d%02d%02d00_NACC2M/'%(today.year,today.month,today.day)
    src_P='http://meteo.caearte.conae.gov.ar/GFS_025/DATA/%d%02d%02d00_ACC/'%(today.year,today.month,today.day)
    os.system('rm '+FORECAST_FOLDER+'* -Rf')
    os.mkdir( FORECAST_TRH_FOLDER)
    os.mkdir( FORECAST_P_FOLDER)
    downloadAll(src_T_RH,FORECAST_TRH_FOLDER)#Temperature and rh
    downloadAll(src_P,FORECAST_P_FOLDER)#precipitation
    renameAll([FORECAST_TRH_FOLDER,FORECAST_P_FOLDER])


def downloadData(start_date,end_date):
    logging.info('Downloading GDAS(fnl)')
    gdas_lib.downloadData(start_date,end_date,GDAS_FOLDER)
    logging.info('Downloading GDAS(anl)')
    gdas_lib.downloadYesterdayAnlData(GDAS_FOLDER)
    logging.info('Downloading IMERG')
    imerg_lib.downloadData(start_date,end_date,IMERG_FOLDER)
    logging.info('Downloading forecast')
    downloadForecast()


#TODO: take into account the utc time. ?
def extractDailyDataFromIMERG(lat,lon,a_date):
    nc_filename=IMERG_FOLDER+imerg_lib.getFilename(a_date)
    grp = nc.Dataset(nc_filename)
    lats = grp.variables['lat'][:]
    lons = grp.variables['lon'][:]
    precipitations=grp.variables['precipitationCal']
    return precipitations[(abs(lats-lat)).argmin(),(abs(lons-lon)).argmin()]

#TODO: take into account the utc time.
def extractDailyDataFromGDAS(lat,lon,a_date,folder,FIELDS,typeOfLevel):
    TIMES=['00','06','12','18']
    epsilon=0.2#TODO:avoid this
    fields_values= dict( (field,[]) for field in FIELDS)
    for a_time in TIMES:
        aux_date=a_date
        if(a_time=='18'): aux_date=a_date-datetime.timedelta(days=1)#utc hack
        grib_filename=folder+gdas_lib.getFilename(aux_date,a_time)
        if(not os.path.isfile(grib_filename)):
            logging.warning('%s not found, but keep going anyways'%grib_filename)
            continue
        grbs=pygrib.open(grib_filename)
        for field in FIELDS:
            grb = grbs.select(name=field,typeOfLevel=typeOfLevel)[0]
            #validate lat,lon
            lats, lons = grb.latlons()
            assert lats.min()<=lat<=lats.max() and lons.min()<=lon<=lons.max()
            #extract the data
            data, lats, lons = grb.data(lat1=lat-epsilon,lat2=lat+epsilon,lon1=lon-epsilon,lon2=lon+epsilon)#TODO:use lat,lon to fabricate lat1,lat2,lon1,lon2
            value=data[0,0]
            if(grb['units']=='K'): value-=273.15 #Kelvin->C
            fields_values[field]+=[ value ]#check this!
    #day ended
    return fields_values

def extractPresentData(lat,lon,start_date,end_date,out_filename):
    output=''
    if(not os.path.isfile(out_filename)): output='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h' + '\n'
    for a_date in daterange(start_date,end_date):
        FIELDS=['2 metre temperature','Relative humidity']
        fields_values=extractDailyDataFromGDAS(lat,lon+360.,a_date,GDAS_FOLDER,FIELDS,typeOfLevel='heightAboveGround')#not sure why it does allow lat to be negative but not lon
        min_T,max_T=np.min(fields_values[FIELDS[0]]),np.max(fields_values[FIELDS[0]])
        mean_T=(min_T+max_T)/2.
        mean_rh=(np.min(fields_values[FIELDS[1]])+np.max(fields_values[FIELDS[1]]))/2.

        precipitation=extractDailyDataFromIMERG(lat,lon,a_date)
        output+=a_date.strftime('%Y-%m-%d')+', '+', '.join([str(min_T),str(mean_T),str(max_T),str(precipitation),str(mean_rh) ]) + ',,'+'\n'
    open(out_filename,'a').write(output)

def extractForecastData(lat,lon,out_filename):
    output='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h' + '\n'
    today=datetime.date.today()
    for a_date in daterange(today,today+datetime.timedelta(hours=168)):
        FIELDS=['2 metre temperature','Relative humidity']
        fields_values=extractDailyDataFromGDAS(lat,lon,a_date,FORECAST_TRH_FOLDER,FIELDS,typeOfLevel='heightAboveGround')#not sure why it does allow lat to be negative but not lon)
        min_T,max_T=np.min(fields_values[FIELDS[0]]),np.max(fields_values[FIELDS[0]])
        mean_T=(min_T+max_T)/2.
        mean_rh=(np.min(fields_values[FIELDS[1]])+np.max(fields_values[FIELDS[1]]))/2.

        fields_values=extractDailyDataFromGDAS(lat,lon,a_date,FORECAST_P_FOLDER,['Total Precipitation'],typeOfLevel='surface')
        precipitation=fields_values['Total Precipitation']
        output+=a_date.strftime('%Y-%m-%d')+', '+', '.join([str(min_T),str(mean_T),str(max_T),str(np.sum(precipitation)),str(mean_rh) ]) + ',,'+'\n'
        print(output)
    open(out_filename.replace('.csv','.forecast.csv'),'w').write(output)

def extractData(params):
    lat,lon,start_date,end_date,out_filename=params
    logging.info('Extracting data to %s'%out_filename)
    extractPresentData(lat,lon,start_date,end_date,out_filename)
    extractForecastData(lat,lon,out_filename)

def joinFullWeather():
    for location in getLocations():
        present_data=open(DATA_FOLDER+location+'.csv','r').read()
        forecast_data=open(DATA_FOLDER+location+'.forecast.csv','r').read()
        open(DATA_FOLDER+location+'.full.csv','w').write(present_data+ '\n'.join(forecast_data.split('\n')[1:]))#remove the header of forecast data

if(__name__ == '__main__'):
    FORMAT='%Y-%m-%d'
    start_date,end_date=None,None
    if(len(sys.argv)>2):
        start_date,end_date= datetime.datetime.strptime(sys.argv[1],FORMAT).date(),datetime.datetime.strptime(sys.argv[2],FORMAT).date()
    elif(len(sys.argv)==1):
        today=datetime.date.today()
        yesterday=today-datetime.timedelta(1)
        start_date,end_date= yesterday,today

    downloadData(start_date,end_date)
    config_parser = ConfigParser()
    config_parser.read('resources/get_weather.cfg')
    params=[]
    for location in config_parser.sections():
        lat=float(config_parser.get(location,'lat'))
        lon=float(config_parser.get(location,'lon'))
        params=params+[[lat,lon,start_date,end_date,DATA_FOLDER+location+'.csv']]

    p = mp.Pool(mp.cpu_count() -1)
    vOpt=p.map(extractData, params)

    joinFullWeather()
