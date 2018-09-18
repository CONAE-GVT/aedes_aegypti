import os
import re
import sys
import urllib
import logging
import datetime
import pygrib
import numpy as np
from os import path
import netCDF4 as nc
import http.cookiejar
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
            new_filename=getFilenameForGDAS(a_date_time.date(),'%02d'%a_date_time.hour,f='00')
            os.rename(dst+'/'+filename,dst+'/'+new_filename)

def downloadAll(url,folder):
    html_list=urllib.request.urlopen(url).read().decode('ISO-8859-1')
    filenames=re.findall(r'.*\<a href=\"([A-Za-z0-9\+\._]+)\"\>.*',html_list)
    for filename in filenames:
        if(filename.endswith('.grib2_NACC_2m') or filename.endswith('.grib2_ACC')):
            open(folder+'/'+filename, 'wb').write(urllib.request.urlopen(url+'/'+filename).read())

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
    downloadDataFromGDAS(start_date,end_date,GDAS_FOLDER)
    logging.info('Downloading IMERG')
    downloadDataFromIMERG(start_date,end_date,IMERG_FOLDER)
    logging.info('Downloading forecast')
    downloadForecast()

#IMERG#
#TODO: take into account the utc time. ?
def getFilenameForIMERG(a_date):
    return '3B-DAY-E.MS.MRG.3IMERG.{year}{month:02}{day:02}-S000000-E235959.V05.nc4'.format(year=a_date.year,month=a_date.month,day=a_date.day)

#https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
def downloadDataFromIMERG(start_date,end_date,folder):
    config_parser = ConfigParser()
    config_parser.read('resources/passwords.cfg')
    passman = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    passman.add_password(None, 'https://urs.earthdata.nasa.gov',config_parser.get('IMERG','username'), config_parser.get('IMERG','password'))
    opener = urllib.request.build_opener(urllib.request.HTTPBasicAuthHandler(passman),urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    urllib.request.install_opener(opener)
    for a_date in daterange(start_date,end_date):
        filename=getFilenameForIMERG(a_date)
        if(path.isfile(folder+'/'+filename)): continue#TODO: Also check filesize
        url='https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGDE.05/{year}/{month:02}/{filename}'.format(year=a_date.year,month=a_date.month,filename=filename)
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request)
        handle = open(folder+'/'+filename, 'wb').write(response.read())

def extractDailyDataFromIMERG(lat,lon,a_date):
    nc_filename=IMERG_FOLDER+getFilenameForIMERG(a_date)
    grp = nc.Dataset(nc_filename)
    lats = grp.variables['lat'][:]
    lons = grp.variables['lon'][:]
    precipitations=grp.variables['precipitationCal']
    return precipitations[(abs(lats-lat)).argmin(),(abs(lons-lon)).argmin()]

#GDAS#
#TODO: take into account the utc time.
def getFilenameForGDAS(a_date,a_time,f):
    return 'gdas1.fnl0p25.{year}{month:02}{day:02}{time}.f{forecast}.grib2'.format(year=a_date.year,month=a_date.month,day=a_date.day,time=a_time,forecast=f)

def downloadDataFromGDAS(start_date,end_date,folder):
    for a_date in daterange(start_date,end_date):
        for a_time in ['00','06','12','18']:
            for a_forecast in ['00','03','06','09']:#a forcast time
                url='http://nomads.ncep.noaa.gov/cgi-bin/filter_gdas_0p25.pl?file=gdas.t%sz.pgrb2.0p25.f0%s&lev_2_m_above_ground=on&var_GUST=on&var_RH=on&var_TCDC=on&var_TMAX=on&var_TMIN=on&var_TMP=on&subregion=&leftlon=-68&rightlon=-60&toplat=-28&bottomlat=-36&dir=%%2Fgdas.%d%02d%02d'%(a_time,a_forecast,a_date.year,a_date.month,a_date.day)
                filename=getFilenameForGDAS(a_date,a_time,f=a_forecast)
                open(folder+'/'+filename, 'wb').write(urllib.request.urlopen(url).read())

def extractDailyDataFromGDAS(lat,lon,a_date,folder,FIELDS,typeOfLevel,f):
    TIMES=['00','06','12','18']
    epsilon=0.5#TODO:avoid this
    fields_values= dict( (field,[]) for field in FIELDS)
    for a_time in TIMES:
        grib_filename=folder+getFilenameForGDAS(a_date,a_time,f)
        if(a_time=='00'): grib_filename=folder+getFilenameForGDAS(a_date + datetime.timedelta(days=1),a_time,f)#utc hack
        if(not os.path.isfile(grib_filename)):
            logging.warning('%s not found, but keep going anyways'%grib_filename)
            continue
        grbs=pygrib.index(grib_filename,'name','typeOfLevel')
        for field in FIELDS:
            grb = grbs.select(name=field,typeOfLevel=typeOfLevel)[0]
            assert (grb.validDate - datetime.timedelta(hours=3,seconds=1)).date() == a_date, '%s vs %s for %s'%( (grb.validDate - datetime.timedelta(hours=3,seconds=1)).date(),a_date,grib_filename) #the second is because I want to take into account the 00 of the next day
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

def extractHistoricData(lat,lon,start_date,end_date,out_filename):
    output=''
    if(not os.path.isfile(out_filename)): output='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h' + '\n'
    for a_date in daterange(start_date,end_date):
        FIELDS=['Minimum temperature','Maximum temperature','2 metre relative humidity']#'2 metre relative humidity' -->'Relative humidity' in the python3 migration
        #to validate that the +360 was ok: 1) gdal_translate a grib to a tif and open qgis with google map as background. 2) use https://www.latlong.net/Show-Latitude-Longitude.html 3)explore.py
        fields_values=extractDailyDataFromGDAS(lat,lon+360.,a_date,GDAS_FOLDER,FIELDS,typeOfLevel='heightAboveGround',f='03')
        min_T,max_T=np.min(fields_values[FIELDS[0]]),np.max(fields_values[FIELDS[1]])
        mean_T=(min_T+max_T)/2.
        mean_rh=(np.min(fields_values[FIELDS[2]])+np.max(fields_values[FIELDS[2]]))/2.

        precipitation=extractDailyDataFromIMERG(lat,lon,a_date)
        output+=a_date.strftime('%Y-%m-%d')+', '+', '.join([str(min_T),str(mean_T),str(max_T),str(precipitation),str(mean_rh) ]) + ',,'+'\n'
    open(out_filename,'a').write(output)

def extractForecastData(lat,lon,out_filename):
    output='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h' + '\n'
    today=datetime.date.today()
    for a_date in daterange(today,today+datetime.timedelta(hours=168)):
        FIELDS=['2 metre temperature','2 metre relative humidity']#'2 metre relative humidity' -->'Relative humidity' in the python3 migration
        fields_values=extractDailyDataFromGDAS(lat,lon,a_date,FORECAST_TRH_FOLDER,FIELDS,typeOfLevel='heightAboveGround',f='00')#not sure why it does allow lat to be negative but not lon)
        min_T,max_T=np.min(fields_values[FIELDS[0]]),np.max(fields_values[FIELDS[0]])
        mean_T=(min_T+max_T)/2.
        mean_rh=(np.min(fields_values[FIELDS[1]])+np.max(fields_values[FIELDS[1]]))/2.

        fields_values=extractDailyDataFromGDAS(lat,lon,a_date,FORECAST_P_FOLDER,['Total Precipitation'],typeOfLevel='surface',f='00')
        precipitation=fields_values['Total Precipitation']
        output+=a_date.strftime('%Y-%m-%d')+', '+', '.join([str(min_T),str(mean_T),str(max_T),str(np.sum(precipitation)),str(mean_rh) ]) + ',,'+'\n'
    open(out_filename.replace('.csv','.forecast.csv'),'w').write(output)

def extractData(params):
    lat,lon,start_date,end_date,out_filename=params
    logging.info('Extracting data to %s'%out_filename)
    #if(os.path.isfile(out_filename)): removeLastLine(out_filename)
    extractHistoricData(lat,lon,start_date,end_date,out_filename)
    extractForecastData(lat,lon,out_filename)

def joinFullWeather():
    for location in getLocations():
        historic_data=open(DATA_FOLDER+location+'.csv','r').read()
        forecast_data=open(DATA_FOLDER+location+'.forecast.csv','r').read()
        open(DATA_FOLDER+location+'.full.csv','w').write(historic_data+ '\n'.join(forecast_data.split('\n')[1:]))#remove the header of forecast data

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
