#based on
import os
import re
import sys
import json
import time
import urllib
import logging
import getpass
import datetime
import http.cookiejar
from utils import daterange
from configparser import ConfigParser

base='https://rda.ucar.edu/apps/'
config_parser = ConfigParser()
config_parser.read('resources/passwords.cfg')
username=config_parser.get('GDAS','username')
password=config_parser.get('GDAS','password')
control_template_filename='resources/ds083.3_control_file.template'
loginurl='https://rda.ucar.edu/cgi-bin/login'
FILENAME_FORMAT='gdas1.fnl0p25.%d%02d%02d%s.f%s.grib2'
LOG_FILENAME='logs/get_weather.log'
logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s',filename=LOG_FILENAME,level=logging.DEBUG)

def init():
    passman = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    passman.add_password(None,'https://rda.ucar.edu', username, password)
    opener = urllib.request.build_opener(urllib.request.HTTPBasicAuthHandler(passman),urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
    urllib.request.install_opener(opener)

def submit(start_date,end_date):
    control_file= open(control_template_filename).read().format(start_date=start_date.strftime('%Y%m%d'+'0000'),end_date=end_date.strftime('%Y%m%d'+'0000')).split('\n')
    control_params={}
    for line in control_file:
        if line.startswith('#') or line=='':
            continue
        li=line.rstrip()
        (key,value)=li.split('=',2)
        control_params[key]=value

    jsondata=json.dumps(control_params).encode('utf-8')
    response=urllib.request.urlopen(urllib.request.Request(base+'request',jsondata,{'Content-type': 'application/json'}) ).read().decode('ISO-8859-1')
    index=re.findall(r'Index[\ ]+:[\ ]+([0-9]+)',response.replace('\n',';'))[0]
    return index

def getStatus(index):
    response = urllib.request.urlopen(base+'/request/'+index).read().decode('ISO-8859-1')
    status=re.findall(r'RequestStatus:[\ ]+([^;]+)',response.replace('\n',';'))[0]
    return status

def waitFor(index):
    for i in range(6,14):#maximum wait is 2**i. (~ 4.5 hr)
        status=getStatus(index)
        if 'Online' in status:
            return
        else:
            logging.info('Waiting %s mins. for %s to be online.'%(2**i /60., index))
            time.sleep(2**i)
    raise GDASError()

def login():
    authdata='email='+username+'&password='+password+'&action=login'
    return urllib.request.urlopen(loginurl,authdata.encode('utf-8')).read()


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
                    open(outpath, 'wb').write( urllib.request.urlopen(remote_filename).read() )
                    logging.info('Download: %s'% remote_filename)

def download(index,folder):
    file_list=json.loads( urllib.request.urlopen(base+'/request/'+index+'/filelist').read() )
    download_files(file_list,folder)

def purge(index):#method not tested
    init()#hack. find a way to avoid this
    request = urllib.request.Request(base+'/request/'+index)
    request.get_method = lambda: 'DELETE'
    logging.info('Purge: %s'% urllib.request.urlopen(request).read())


def getFilename(a_date,a_time,f):
    return FILENAME_FORMAT%(a_date.year,a_date.month,a_date.day,a_time,f)

def downloadData(start_date,end_date,folder):
    init()
    index=submit(start_date,end_date)
    waitFor(index)
    download(index,folder)
    purge(index)

class GDASError(Exception):
    pass

if __name__=='__main__':
    FORMAT='%Y-%m-%d'
    start_date,end_date=None,None
    if(len(sys.argv)>2):
        start_date,end_date= datetime.datetime.strptime(sys.argv[1],FORMAT).date(),datetime.datetime.strptime(sys.argv[2],FORMAT).date()
    elif(len(sys.argv)==1):
        start_date,end_date=datetime.date.today()-datetime.timedelta(days=1),datetime.date.today()

    downloadData(start_date,end_date,'data/public/gdas')
