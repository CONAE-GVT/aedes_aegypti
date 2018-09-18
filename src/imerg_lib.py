from os import path
import urllib
import datetime
import http.cookiejar
from utils import daterange
#TODO:put this in a config file
base='https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L3/GPM_3IMERGDE.05/'
username='***REMOVED***'
password='***REMOVED***'

def getFilename(a_date):
    return '3B-DAY-E.MS.MRG.3IMERG.{year}{month:02}{day:02}-S000000-E235959.V05.nc4'.format(year=a_date.year,month=a_date.month,day=a_date.day)

#https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+Python
def downloadData(start_date,end_date,folder):
    passman = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    passman.add_password(None, 'https://urs.earthdata.nasa.gov', username, password)
    opener = urllib.request.build_opener(urllib.request.HTTPBasicAuthHandler(passman),urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    urllib.request.install_opener(opener)
    for a_date in daterange(start_date,end_date):
        filename=getFilename(a_date)
        if(path.isfile(folder+'/'+filename)): continue#TODO: Also check filesize
        url=base+'{year}/{month:02}/{filename}'.format(year=a_date.year,month=a_date.month,filename=filename)
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request)
        handle = open(folder+'/'+filename, 'wb').write(response.read())


if __name__=='__main__':
    downloadData(datetime.date.today()-datetime.timedelta(days=1),datetime.date.today())
