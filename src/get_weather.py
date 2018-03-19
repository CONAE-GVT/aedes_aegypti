from HTMLParser import HTMLParser
from dateutil import relativedelta
import datetime
import httplib2
import time
import sys
import os


url="https://www.wunderground.com/history/airport/SACO/%s/%s/01/MonthlyHistory.html"
OUT_FILENAME='wunderground.csv'

class TableParser(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.inside_table = None
        self.tag=None
        self.csv=''

    def handle_starttag(self, tag, attrs):
        if('table'==tag and ('id', 'obsTable') in attrs):
            self.inside_table=True

        if(self.inside_table):
            self.tag=tag
            for colspan in [int(pair[1]) for pair in attrs if tag in ('th','td') and pair[0]=='colspan']:
                self.csv+=','*(colspan-1)#the -1 is because we are always printing a ',' at the end of th or td

    def handle_endtag(self, tag):
        if self.inside_table:
            if(tag in ('th','td')):
                self.csv+=','
            elif(tag=='tr'):
                self.csv+='\n'

            if('table'==tag):
                self.inside_table=False
            self.tag=None

    def handle_data(self, data):
        if self.inside_table and self.tag in ('th','td','span','a'):
            self.csv+=data.replace(',','/').replace('\n','')



def daterange(start_date, end_date):
    delta=relativedelta.relativedelta(end_date,start_date)
    months= delta.years * 12 + delta.months
    for n in range(months):
        yield start_date + relativedelta.relativedelta(months=n)

def retrieve_data(year,month):
    data_url= url % (year,month)
    headers = {
        "Content-type": "text/plain",
        "Accept": "application/xml"
    }
    headers, response = http.request(data_url, "GET")
    parser = TableParser()
    parser.feed(response)
    content=parser.csv
    content=('\n').join(['%i-%i-'%(year,month)+line if line.split(',')[0].isdigit() else line for line in content.split('\n')])#we add the year and month at the begining

    if(not os.path.isfile(OUT_FILENAME)):
        file = open(OUT_FILENAME,'w')
    else:
        file = open(OUT_FILENAME,'a')
        content='\n'.join(content.split('\n')[2:])

    file.write(content)
    file.close()

#extract just the data we need and put it in an specific order so it can be used by other script.
def prune_data():
    file = open(OUT_FILENAME,'r')
    output='Date,Minimum Temp (C),Mean Temperature (C),Maximum Temp (C),Rain (mm),Relative Humidity %,CloudCover,Mean Wind SpeedKm/h' + '\n'
    for line_number,line in enumerate(file):
        if(line_number<2): continue#skip headers
        fields=line.split(',')
        output+=','.join([fields[0],fields[3],fields[2],     fields[1],         fields[19], fields[8],            '',         fields[17]]) + '\n'

    open(OUT_FILENAME.replace('.csv','_pruned.csv'),'w').write(output)

#main
http = httplib2.Http()
if(__name__ == '__main__'):
    #retrieve_data('2017','03')
    for d in daterange(datetime.date(2017, 7,01),datetime.date(2017, 8,01)):
        #print(str(d.year) + ' '+str(d.month))
        time.sleep(15)
        retrieve_data(d.year,d.month)

    #put just the info we need into a file, that will be ready to be used by otero_precipitation script.
    prune_data()
