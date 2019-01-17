from html.parser import HTMLParser
from dateutil import relativedelta
import datetime
import urllib.request
import time
import sys
import os
import json



OUT_FILENAME='wunderground.csv'
#https://api-ak.wunderground.com/api/d8585d80376a429e/history_2014070120140731/lang:EN/units:metric/bestfct:1/v:2.0/q/SACO.json?showObs=0&ttl=120
url='https://api-ak.wunderground.com/api/d8585d80376a429e/history_{year}{month:02}01{year}{month:02}31/lang:EN/units:metric/bestfct:1/v:2.0/q/SACO.json?showObs=0&ttl=120'




def daterange(start_date, end_date):
    delta=relativedelta.relativedelta(end_date,start_date)
    months= delta.years * 12 + delta.months
    for n in range(months):
        yield start_date + relativedelta.relativedelta(months=n)

#main
if(__name__ == '__main__'):
    output='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h' + '\n'
    open(OUT_FILENAME,'a').write(output)
    for d in daterange(datetime.date(2018, 7,1),datetime.date(2019, 1,10)):
        #print(str(d.year) + ' '+str(d.month))
        json_content=urllib.request.urlopen(url.format(year=d.year,month=d.month)).read()
        #json_content=open('content.json','r').read()
        content=json.loads(json_content)
        for day_info in content['history']['days']:
            summary=day_info['summary']
            year=summary['date']['year']
            month=summary['date']['month']
            day=summary['date']['day']
            min_temperature=summary['min_temperature']
            mean_temperature=summary['temperature']
            max_temperature=summary['max_temperature']
            rain=summary['precip']
            relative_humidity=(summary['min_humidity']+summary['max_humidity'])/2.
            #Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h
            output='%s-%s-%s, %s,%s,%s,%s,%s,,,\n'%(year,month,day,min_temperature,mean_temperature,max_temperature,rain,relative_humidity)
            open(OUT_FILENAME,'a').write(output)

        time.sleep(15)
