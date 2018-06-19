import datetime
import imerg_lib
import gdas_lib
DATA_FOLDER='data/public/'
IMERG_FOLDER=DATA_FOLDER+'/imerg/'
GDAS_FOLDER=DATA_FOLDER+'/gdas/'
OUT_FILENAME=DATA_FOLDER+'weather.csv'

def downloadData(start_date,end_date):
    gdas_lib.downloadData(start_date,end_date,GDAS_FOLDER)
    imerg_lib.downloadData(start_date,end_date,IMERG_FOLDER)

#TODO: take into account the utc time. ?
def extractDailyDataFromIMERG(lat,lon,the_date):
    nc_filename=IMERG_FOLDER+'3B-DAY-E.MS.MRG.3IMERG.%d%02d%02d-S000000-E235959.V05.nc4.nc'%(the_date.year,the_date.month,the_date.day)
    grp = nc.Dataset(nc_filename)
    #print( grp.variables['lat'])
    lats = grp.variables['lat'][:]
    lons = grp.variables['lon'][:]
    precipitations=grp.variables['precipitationCal']
    return precipitations[(abs(lats-lat)).argmin(),(abs(lons-lon)).argmin()]

#TODO: take into account the utc time.
def extractDailyDataFromGDAS(lat,lon,the_date):
    TIMES=['00','06','12','18']
    #FIELDS=['Minimum temperature','Maximum temperature','Relative humidity']
    FIELDS=['Temperature','Relative humidity']
    lon+=360.#not sure why it does allow lat to be negative
    epsilon=0.1#TODO:avoid this
    #print('%s, %s'%(lat+epsilon,lon+epsilon))
    for a_time in TIMES:
        #gdas1.fnl0p25.2017070106.f06.grib2.spasub.aguirre296700
        #grib_filename=folder+'gdas1.fnl0p25.%d%02d%02d%s.f03.grib2.spasub.aguirre296700'%(the_date.year,the_date.month,the_date.day,a_time)
        aux_date=the_date
        #if(a_time=='18'): aux_date=the_date-datetime.timedelta(days=1)#utc hack
        grib_filename=GDAS_FOLDER+'gdas1.fnl0p25.%d%02d%02d%s.f09.grib2.spasub.aguirre298079'%(aux_date.year,aux_date.month,aux_date.day,a_time)
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

def extractData(lat,lon,start_date,end_date):
    file = open(OUT_FILENAME,'w')
    output='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h' + '\n'
    for a_date in daterange(start_date,end_date):
        min_T,mean_T,max_T,mean_rh=extractDailyDataFromGDAS(lat,lon,a_date)
        rain=extractDailyDataFromGPM(lat,lon,a_date)
        output+=a_date.strftime('%Y-%m-%d')+', '+', '.join([str(min_T),str(mean_T),str(max_T),str(rain),str(mean_rh) ]) + ',,'+'\n'
        print(output)
    open(OUT_FILENAME,'w').write(output)

if(__name__ == '__main__'):
    FORMAT='%Y-%m-%d'
    start_date,end_date= datetime.datetime.strptime(sys.argv[1],FORMAT).date(),datetime.datetime.strptime(sys.argv[2],FORMAT).date()
    downloadData(start_date,end_date)
    lat,lon=-31.420083,-64.188776
    extractData(lat,lon,start_date,end_date)
