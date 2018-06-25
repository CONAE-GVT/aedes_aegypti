import sys
import logging
import datetime
import gdas_lib
import imerg_lib
import pygrib
import numpy as np
import netCDF4 as nc
from utils import daterange

DATA_FOLDER='data/public/'
IMERG_FOLDER=DATA_FOLDER+'/imerg/'
GDAS_FOLDER=DATA_FOLDER+'/gdas/'
OUT_FILENAME=DATA_FOLDER+'weather.csv'
LOG_FILENAME='logs/get_weather.log'
logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s',filename=LOG_FILENAME,level=logging.DEBUG)

def downloadData(start_date,end_date):
    logging.info('Downloading GDAS')
    gdas_lib.downloadData(start_date,end_date,GDAS_FOLDER)
    logging.info('Downloading IMERG')
    imerg_lib.downloadData(start_date,end_date,IMERG_FOLDER)

#TODO: take into account the utc time. ?
def extractDailyDataFromIMERG(lat,lon,a_date):
    nc_filename=IMERG_FOLDER+imerg_lib.getFilename(a_date)
    grp = nc.Dataset(nc_filename)
    lats = grp.variables['lat'][:]
    lons = grp.variables['lon'][:]
    precipitations=grp.variables['precipitationCal']
    return precipitations[(abs(lats-lat)).argmin(),(abs(lons-lon)).argmin()]

#TODO: take into account the utc time.
def extractDailyDataFromGDAS(lat,lon,a_date):
    TIMES=['00','06','12','18']
    FIELDS=['2 metre temperature','Relative humidity']
    lon+=360.#not sure why it does allow lat to be negative
    epsilon=0.1#TODO:avoid this
    fields_values= dict( (field,[]) for field in FIELDS)
    for a_time in TIMES:
        aux_date=a_date
        if(a_time=='18'): aux_date=a_date-datetime.timedelta(days=1)#utc hack
        grib_filename=GDAS_FOLDER+gdas_lib.getFilename(aux_date,a_time)
        grbs=pygrib.open(grib_filename)
        for field in FIELDS:
            grb = grbs.select(name=field,typeOfLevel='heightAboveGround')[0]
            #validate lat,lon
            lats, lons = grb.latlons()
            assert lats.min()<=lat<=lats.max() and lons.min()<=lon<=lons.max()
            #extract the data
            data, lats, lons = grb.data(lat1=lat-epsilon,lat2=lat+epsilon,lon1=lon-epsilon,lon2=lon+epsilon)#TODO:use lat,lon to fabricate lat1,lat2,lon1,lon2
            value=data[0,0]
            if(grb['units']=='K'): value-=273.15 #Kelvin->C
            fields_values[field]+=[ value ]#check this!
    #day ended
    min_T=np.min(fields_values[FIELDS[0]])
    max_T=np.max(fields_values[FIELDS[0]])
    mean_T=(min_T+max_T)/2.
    mean_rh=(np.min(fields_values[FIELDS[1]])+np.max(fields_values[FIELDS[1]]))/2.
    return min_T,mean_T,max_T,mean_rh

def extractData(lat,lon,start_date,end_date):
    output='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h' + '\n'
    for a_date in daterange(start_date,end_date):
        min_T,mean_T,max_T,mean_rh=extractDailyDataFromGDAS(lat,lon,a_date)
        rain=extractDailyDataFromIMERG(lat,lon,a_date)
        output+=a_date.strftime('%Y-%m-%d')+', '+', '.join([str(min_T),str(mean_T),str(max_T),str(rain),str(mean_rh) ]) + ',,'+'\n'
    open(OUT_FILENAME,'w').write(output)

if(__name__ == '__main__'):
    FORMAT='%Y-%m-%d'
    start_date,end_date= datetime.datetime.strptime(sys.argv[1],FORMAT).date(),datetime.datetime.strptime(sys.argv[2],FORMAT).date()
    downloadData(start_date,end_date)
    lat,lon=-31.420083,-64.188776
    extractData(lat,lon,start_date,end_date)
