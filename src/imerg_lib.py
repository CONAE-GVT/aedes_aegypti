import urllib2
import datetime
import cookielib

base='https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGDE.05/'
username='***REMOVED***'
password='***REMOVED***'

#different version than in get_weather
def daterange(start_date, end_date):
    current_date=start_date
    while current_date<=end_date:
        yield current_date
        current_date=current_date+datetime.timedelta(days=1)


#https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
def downloadData(start_date,end_date,folder):
    passman = urllib2.HTTPPasswordMgrWithDefaultRealm()
    passman.add_password(None, 'https://urs.earthdata.nasa.gov', username, password)
    opener = urllib2.build_opener(urllib2.HTTPBasicAuthHandler(passman),urllib2.HTTPCookieProcessor(cookielib.CookieJar()))
    urllib2.install_opener(opener)
    for a_date in daterange(start_date,end_date):
        filename='3B-DAY-E.MS.MRG.3IMERG.{year}{month:02}{day:02}-S000000-E235959.V05.nc4'.format(year=a_date.year,month=a_date.month,day=a_date.day)
        url=base+'{year}/{month:02}/{filename}'.format(year=a_date.year,month=a_date.month,filename=filename)
        request = urllib2.Request(url)
        response = urllib2.urlopen(request)
        handle = open(folder+'/'+filename, 'w').write(response.read())


if __name__=='__main__':
    downloadData(datetime.date.today()-datetime.timedelta(days=1),datetime.date.today())
