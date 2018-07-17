#https://jswhit.github.io/pygrib/docs/
import os
import sys
import pygrib
import datetime
import numpy as np
import netCDF4 as nc
from configparser import ConfigParser

DATA_FOLDER='data/public/'
OUT_FILENAME=DATA_FOLDER+'weather.csv'

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + datetime.timedelta(n)

def extractDataWgrib2(lat,lon,grib_filename,epsilon):
    os.system('wgrib2 {grib_filename} -undefine out-box {lon_min}:{lon_max} {lat_min}:{lat_max} -csv tmp.csv > /dev/null'.format(grib_filename=grib_filename,lon_min=lon-epsilon,lon_max=lon+epsilon,lat_min=lat-epsilon,lat_max=lat+epsilon))
    fields_values=dict()
    for line in open('tmp.csv','r').readlines():
        values=line.split(',')
        fields_values[values[2].replace('"','')]=float(values[-1])
    return fields_values

def extractDailyDataFromReanalysis(lat,lon,grib_filename):
    FIELDS=['Minimum temperature','2 metre temperature','Maximum temperature','Relative humidity']
    lon+=360.#not sure why it does allow lat to be negative
    epsilon=0.1#TODO:avoid this
    grbs=pygrib.open(grib_filename)
    fields_values= dict( (field,[]) for field in FIELDS)
    for field in FIELDS:
        grb = grbs.select(name=field,typeOfLevel='heightAboveGround')[0]
        #validate lat,lon
        lats, lons = grb.latlons()
        assert lats.min()<=lat<=lats.max() and lons.min()<=lon<=lons.max()
        #extract the data
        data, lats, lons = grb.data(lat1=lat-epsilon,lat2=lat+epsilon,lon1=lon-epsilon,lon2=lon+epsilon)#TODO:use lat,lon to fabricate lat1,lat2,lon1,lon2
        #if(not data): print(extractDataWgrib2(lat,lon-360.,grib_filename,epsilon))
        fields_values[field]=data[0,0]
    #day ended
    wgrib2_field_values=extractDataWgrib2(lat,lon-360.,grib_filename,epsilon)
    #FIELDS=['Minimum temperature','2 metre temperature','Maximum temperature','Relative humidity']
    assert abs(fields_values['Minimum temperature']-wgrib2_field_values['TMIN'])<1e-3,'TMIN'
    assert abs(fields_values['2 metre temperature']-wgrib2_field_values['TMP'])<1e-3,'TMP'
    assert abs(fields_values['Maximum temperature']-wgrib2_field_values['TMAX'])<1e-3,'TMAX'
    assert abs(fields_values['Relative humidity']-wgrib2_field_values['RH'])<1e-3,'RH'

    min_T=np.min(fields_values[FIELDS[0]])
    max_T=np.max(fields_values[FIELDS[0]])
    mean_T=(min_T+max_T)/2.
    mean_rh=(np.min(fields_values[FIELDS[1]])+np.max(fields_values[FIELDS[1]]))/2.
    return min_T,mean_T,max_T,mean_rh


if(__name__ == '__main__'):
    folder='data/public/gdas/'
    config_parser = ConfigParser()
    config_parser.read('resources/get_weather.cfg')
    for location in config_parser.sections():
        if(location in ['rio_cuarto','marcos_juarez','la_para','salsipuedes']): continue
        print('location: %s...'%location),
        lat=float(config_parser.get(location,'lat'))
        lon=float(config_parser.get(location,'lon'))
        for filename in os.listdir(folder):
            if (filename.endswith('f03.grib2') ):
                extractDailyDataFromReanalysis(lat,lon,folder+filename)
        print('ok')
