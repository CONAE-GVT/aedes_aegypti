#based on
###################################################################################
###
###     Title: rdams-client.py
###    Author: Doug Schuster, schuster@ucar.edu
###      Date: 10/24/2013
###   Purpose: List dataset metadata, subset data subset requests, check on request
###            status.
###
###  SVN File: $HeadURL: https://subversion.ucar.edu/svndss/schuster/rest_client/rdams-client.py $
####################################################################################

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
cookie_file='auth.rda_ucar_edu'
loginurl='https://rda.ucar.edu/cgi-bin/login'
FILENAME_FORMAT='gdas1.fnl0p25.%d%02d%02d%s.f09.grib2'
LOG_FILENAME='logs/get_weather.log'
logging.basicConfig(format='%(levelname)s: %(asctime)s %(message)s',filename=LOG_FILENAME,level=logging.DEBUG)

def getOpenedRequest(url,data=None,data_type=None):
    opener = add_http_auth(url,username,password)
    request = None
    if(data is not None and data_type is not None):
        request=urllib2.Request(url,data,data_type)
    else:
        request=urllib2.Request(url)
    opened_request=None
    try:
    	opened_request = opener.open(request)
    except urllib2.HTTPError, e:
    	if e.code == 401:
    		print 'RDA username and password invalid.\n'
    		sys.exit()
    return opened_request

# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update_progress(progress,outdir):
    barLength = 20 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n\n"
    block = int(round(barLength*progress))
    text = "\rDownloading Request to './{0}' directory.  Download Progress: [{1}] {2}% {3}".format( outdir,"="*block + " "*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

# download_file(remfile,outfile) : download a file from a remote server (remfile) to a local location (outfile)
def download_file(remfile,outfile):
    frequest = urllib2.Request(remfile)
    fresponse = urllib2.urlopen(remfile)
    handle = open(outfile, 'w')
    handle.write(fresponse.read())
    handle.close()

# add_http_auth(url,user,pasw): add authentication information to opener and return opener
def add_http_auth(url,user,pasw):
        passman = urllib2.HTTPPasswordMgrWithDefaultRealm()
        passman.add_password(None, url, username, password)
        authhandler = urllib2.HTTPBasicAuthHandler(passman)
        opener = urllib2.build_opener(authhandler)
        urllib2.install_opener(opener)
        return opener

# add_http_cookie(url,authstring): Get and add authentication cookie to http file download handler
def add_http_cookie(url,authstring):
        cj = cookielib.MozillaCookieJar(cookie_file)
        openrf=urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
        frequest = urllib2.Request(url,authstring)
        cj.add_cookie_header(frequest)
        response=openrf.open(frequest)
        openerf = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
        urllib2.install_opener(openerf)

# download_files(filelist,directory): Download multiple files from the rda server and save them to a local directory
def download_files(filelist,directory):
        backslash='/'
        filecount=0
        percentcomplete=0
        localsize=''
        length=0
        length=len(filelist)
        if not os.path.exists(directory):
                os.makedirs(directory)
        for key, value in filelist.iteritems():
                downloadpath,localfile=key.rsplit("/",1)
                localfile='.'.join(localfile.split('.')[:-2])#gdas1.fnl0p25.2017090212.f03.grib2.spasub.aguirre296700-->gdas1.fnl0p25.2017090212.f03.grib2
                outpath=directory+backslash+localfile
                percentcomplete=(float(filecount)/float(length))
                update_progress(percentcomplete,directory)
                if os.path.isfile(outpath):
                        localsize=os.path.getsize(outpath)
                        if(str(localsize) != value):
                                download_file(key,outpath)
                elif(not os.path.isfile(outpath)):
                        download_file(key,outpath)

                filecount=filecount+1
                percentcomplete=(float(filecount)/float(length))
        update_progress(percentcomplete,directory)

def getFilename(a_date,a_time):
    return 'gdas1.fnl0p25.%d%02d%02d%s.f09.grib2'%(a_date.year,a_date.month,a_date.day,a_time)

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
    opened_request=getOpenedRequest(base+'request',jsondata,{'Content-type': 'application/json'})
    response=opened_request.read()
    index=re.findall(r'Index[\ ]+:[\ ]+([0-9]+)',response.replace('\n',';'))[0]
    return index

def getStatus(index):
    response=getOpenedRequest(base+'/request/'+index).read()
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

def download(index,folder):
    # get cookie required to download data files
    file_list=json.loads( getOpenedRequest(base+'/request/'+index+'/filelist').read() )

    authdata='email='+username+'&password='+password+'&action=login'
    add_http_cookie(loginurl,authdata)#this change the opener, so it needs to be after the last getOpenedRequest (which also changes the opener)
    download_files(file_list,folder)

def purge(index):
    url=base+'/request/'+index
    opener = add_http_auth(url,username,password)
    request = urllib2.Request(url)
    request.get_method = lambda: 'DELETE'
    opened_request = opener.open(request)
    print(opened_request.read())

def downloadData(start_date,end_date,folder):
    #just download the data we don't already have
    for a_date in daterange(start_date,end_date):
        if(os.path.isfile(folder+getFilename(a_date,'00'))): start_date=a_date

    index=submit(start_date,end_date)
    waitFor(index)
    download(index,folder)
    purge(index)

if __name__=='__main__':
    downloadData(datetime.date.today()-datetime.timedelta(days=1),datetime.date.today())
