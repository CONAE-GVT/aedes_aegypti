#https://jswhit.github.io/pygrib/docs/
import numpy as np
import datetime
import pygrib
import netCDF4 as nc
import sys

DATA_FOLDER='data/public/'
OUT_FILENAME=DATA_FOLDER+'weather.csv'

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + datetime.timedelta(n)

#TODO: take into account the utc time. ?
def extractDailyDataFromGPM(lat,lon,the_date):
    folder=DATA_FOLDER+'gpm/'
    nc_filename=folder+'3B-DAY-E.MS.MRG.3IMERG.%d%02d%02d-S000000-E235959.V05.nc4.nc'%(the_date.year,the_date.month,the_date.day)
    grp = nc.Dataset(nc_filename)
    #print( grp.variables['lat'])
    lats = grp.variables['lat'][:]
    lons = grp.variables['lon'][:]
    precipitations=grp.variables['precipitationCal']
    return precipitations[(abs(lats-lat)).argmin(),(abs(lons-lon)).argmin()]

#TODO: take into account the utc time.

def extractDailyDataFromReanalysis(lat,lon,the_date):
    TIMES=['00','06','12','18']
    #FIELDS=['Minimum temperature','Maximum temperature','Relative humidity']
    FIELDS=['Temperature','Relative humidity']
    folder=DATA_FOLDER+'reanalysis/'
    lon+=360.#not sure why it does allow lat to be negative
    epsilon=0.1#TODO:avoid this
    #print('%s, %s'%(lat+epsilon,lon+epsilon))
    for a_time in TIMES:
        #gdas1.fnl0p25.2017070106.f06.grib2.spasub.aguirre296700
        #grib_filename=folder+'gdas1.fnl0p25.%d%02d%02d%s.f03.grib2.spasub.aguirre296700'%(the_date.year,the_date.month,the_date.day,a_time)
        aux_date=the_date
        #if(a_time=='18'): aux_date=the_date-datetime.timedelta(days=1)#utc hack


        grib_filename=folder+'gdas1.fnl0p25.%d%02d%02d%s.f09.grib2.spasub.aguirre298079'%(aux_date.year,aux_date.month,aux_date.day,a_time)
        grbs=pygrib.open(grib_filename)
        fields_values= dict( (field,[]) for field in FIELDS)
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

'''
def extractDailyDataFromReanalysis(lat,lon,the_date):
    folder=DATA_FOLDER+'reanalysis/'
    TIMES=['00','06','12','18']
    FIELDS=['TMP_P0_L1_GLL0','RH_P0_L103_GLL0']
    fields_values= dict( (field,[]) for field in FIELDS)
    for a_time in TIMES:
        aux_date=the_date
        nc_filename=folder+'gdas1.fnl0p25.%d%02d%02d%s.f09.grib2.aguirre298213.nc'%(aux_date.year,aux_date.month,aux_date.day,a_time)
        grp = nc.Dataset(nc_filename)
        #print(grp.variables)
        lats = grp.variables['lat_0'][:]
        lons = grp.variables['lon_0'][:]
        lat_index=(abs(lats-lat)).argmin()
        lon_index=(abs(lons-lon)).argmin()

        #print(grp.variables['R_H_L107'][:,0,(,])
        #print(grp.variables['time'][:])
        #min_T=0
        #max_T=0
        #mean_T=np.mean(grp.variables['TMP_L107'][:,2,lat_index,lon_index]-273.15)
        #mean_rh=np.mean(grp.variables['R_H_L107'][:,2,lat_index,lon_index])

        for field in FIELDS:
            value=grp.variables[field][:,lat_index,lon_index]
            if(grp.variables[field].units=='K'): value-=273.15
            fields_values[field]+=[ value ]#check this!
    #print(fields_values)
    min_T=0
    max_T=0
    mean_T=np.mean(fields_values['TMP_P0_L1_GLL0'])
    mean_rh=np.mean(fields_values['RH_P0_L103_GLL0'])
    return min_T,mean_T,max_T,mean_rh
'''
def extractData(lat,lon,start_date,end_date):
    file = open(OUT_FILENAME,'w')
    output='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h' + '\n'
    for a_date in daterange(start_date,end_date):
        min_T,mean_T,max_T,mean_rh=extractDailyDataFromReanalysis(lat,lon,a_date)
        rain=extractDailyDataFromGPM(lat,lon,a_date)
        output+=a_date.strftime('%Y-%m-%d')+', '+', '.join([str(min_T),str(mean_T),str(max_T),str(rain),str(mean_rh) ]) + ',,'+'\n'
        print(output)
    open(OUT_FILENAME,'w').write(output)

###just to test###
import matplotlib.pyplot as plt
def show(the_date):
    #precipitations
    folder=DATA_FOLDER+'gpm/'
    nc_filename=folder+'3B-DAY-E.MS.MRG.3IMERG.%d%02d%02d-S000000-E235959.V05.nc4.nc'%(the_date.year,the_date.month,the_date.day)
    grp = nc.Dataset(nc_filename)
    data=grp.variables['precipitationCal']
    plt.imshow(data)
    #plt.show()

    #Temperature
    folder=DATA_FOLDER+'reanalysis/'
    grib_filename=folder+'gdas1.fnl0p25.%d%02d%02d%s.f09.grib2.spasub.aguirre296700'%(the_date.year,the_date.month,the_date.day,'18')
    grbs=pygrib.open(grib_filename)
    grb = grbs.select(name='Maximum temperature',typeOfLevel='heightAboveGround')[0]
    lats, lons = grb.latlons()
    data, lats, lons = grb.data(lat1=lats.min(),lat2=lats.max(),lon1=lons.min(),lon2=lons.max())#TODO:use lat,lon to fabricate lat1,lat2,lon1,lon2
    plt.figure()
    plt.imshow(data)
    plt.show()
###just to test###
import utils
def compare(start_date,end_date):
    start_date=datetime.date(2017, 7, 2)
    end_date=datetime.date(2018, 5, 1)
    temperatures_saco= np.array(utils.getAverageTemperaturesFromCsv('data/public/wunderground_SACO.csv',start_date,end_date))
    temperatures_gdas= np.array(utils.getAverageTemperaturesFromCsv('data/public/weather.csv',start_date,end_date))
    plt.plot(temperatures_saco,'-m', label='SACO')
    plt.plot(temperatures_gdas,'-r', label='gdas')

    delta_temperature=temperatures_saco-temperatures_gdas
    plt.plot(delta_temperature,'-b', label='SACO-gdas')
    print('%s +/- %s'%(np.mean(delta_temperature),np.std(delta_temperature)) )

    delta_temperature=temperatures_saco-temperatures_gdas
    plt.plot(delta_temperature,'-b', label='SACO-gdas')
    print('%s +/- %s'%(np.mean(delta_temperature),np.std(delta_temperature)) )
    plt.legend(loc=0)
    plt.show()



if(__name__ == '__main__'):
    FORMAT='%Y-%m-%d'
    start_date,end_date= datetime.datetime.strptime(sys.argv[1],FORMAT).date(),datetime.datetime.strptime(sys.argv[2],FORMAT).date()
    lat,lon=-31.420083,-64.188776
    extractData(lat,lon,start_date,end_date)
    #show(start_date)
    compare(start_date,end_date)
