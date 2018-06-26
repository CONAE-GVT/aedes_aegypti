#based on
import os
import re
import sys
import json
import time
import logging
import urllib2
import getpass
import datetime
import cookielib
from utils import daterange

base='https://rda.ucar.edu/apps/'
username='***REMOVED***'
password='***REMOVED***'
control_template_filename='resources/ds083.3_control_file.template'
loginurl='https://rda.ucar.edu/cgi-bin/login'
FILENAME_FORMAT='gdas1.fnl0p25.%d%02d%02d%s.f09.grib2'
LOG_FILENAME='logs/get_weather.log'
logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s',filename=LOG_FILENAME,level=logging.DEBUG)

def init():
    passman = urllib2.HTTPPasswordMgrWithDefaultRealm()
    passman.add_password(None,'https://rda.ucar.edu', username, password)
    opener = urllib2.build_opener(urllib2.HTTPBasicAuthHandler(passman),urllib2.HTTPCookieProcessor(cookielib.CookieJar()))
    urllib2.install_opener(opener)

def submit(start_date,end_date):
    control_file= open(control_template_filename).read().format(start_date=start_date.strftime('%Y%m%d'+'0000'),end_date=end_date.strftime('%Y%m%d'+'0000')).split('\n')
    control_params={}
    for line in control_file:
        if line.startswith('#') or line=='':
            continue
        li=line.rstrip()
        (key,value)=li.split('=',2)
        control_params[key]=value
    jsondata='{'
    for k in control_params.keys():
        jsondata+='"'+k+'"'+":"+'"'+control_params[k]+'",'
    jsondata = jsondata[:-1]
    jsondata+='}'

    response=urllib2.urlopen(urllib2.Request(base+'request',jsondata,{'Content-type': 'application/json'}) ).read()
    index=re.findall(r'Index[\ ]+:[\ ]+([0-9]+)',response.replace('\n',';'))[0]
    return index

def getStatus(index):
    response = urllib2.urlopen(base+'/request/'+index).read()
    status=re.findall(r'RequestStatus:[\ ]+([^;]+)',response.replace('\n',';'))[0]
    return status

def waitFor(index):
    for i in range(6,15):#maximum wait is 2**i. (~ 9 hr)
        status=getStatus(index)
        if 'Online' in status:
            break
        else:
            logging.info('Waiting %s mins. for %s to be online.'%(2**i /60., index))
            time.sleep(2**i)

def login():
    authdata='email='+username+'&password='+password+'&action=login'
    return urllib2.urlopen(loginurl,authdata).read()


def download_files(filelist,directory):
        login()
        localsize=''
        if not os.path.exists(directory):
            os.makedirs(directory)
        for remote_filename, remote_filesize in filelist.iteritems():
                downloadpath,local_filename=remote_filename.rsplit("/",1)
                #gdas1.fnl0p25.2017090212.f03.grib2.spasub.aguirre296700-->gdas1.fnl0p25.2017090212.f03.grib2
                local_filename='.'.join(local_filename.split('.')[:-2])
                outpath=directory+'/'+local_filename
                #if the file do not exist or the sizes do not match, download
                is_file_inexistant_or_incomplete=not os.path.isfile(outpath) or (os.path.isfile(outpath) and str(os.path.getsize(outpath)) !=remote_filesize)
                if is_file_inexistant_or_incomplete:
                    #downloadthe file
                    open(outpath, 'w').write( urllib2.urlopen(remote_filename).read() )
                    logging.info('Download: %s'% remote_filename)

def download(index,folder):
    file_list=json.loads( urllib2.urlopen(base+'/request/'+index+'/filelist').read() )
    download_files(file_list,folder)

def purge(index):
    init()#hack. find a way to avoid this
    request = urllib2.Request(base+'/request/'+index)
    request.get_method = lambda: 'DELETE'
    logging.info('Purge: %s'% urllib2.urlopen(request).read())


def getFilename(a_date,a_time):
    return FILENAME_FORMAT%(a_date.year,a_date.month,a_date.day,a_time)

def downloadData(start_date,end_date,folder):
    init()
    index=submit(start_date,end_date)
    waitFor(index)
    download(index,folder)
    purge(index)

if __name__=='__main__':
    downloadData(datetime.date.today()-datetime.timedelta(days=1),datetime.date.today())
