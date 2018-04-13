import collections
import numpy as np
import pylab as pl
import datetime
import math
import xlrd



#print(getValuesFromXls('data/Indices aedicos Historicos.xlsx',date(2015, 1, 1),0,2,'SE',1))
def getValuesFromXls(filename,initial_date,date_column,value_column,filter_value=None,filter_column=None):#Maybe add a filter
    #Loads only current sheets to memory
    workbook = xlrd.open_workbook(filename, on_demand = True)

    # Load a specific sheet by index
    worksheet = workbook.sheet_by_index(0)

    # Retrieve the data pairs(key_column,value_column)
    data=collections.OrderedDict()
    last_date_value=None
    for row in range(2,worksheet.nrows):

        if(not worksheet.cell(row, value_column).value):
            continue
        if(worksheet.cell(row, date_column).value):
            last_date_value=worksheet.cell(row, date_column).value
        if(not last_date_value):
            continue
        if(filter_value and filter_column and not worksheet.cell(row, filter_column).value==filter_value):
            continue

        #get the days passed from initial_date to the date in the excel
        year,month,day,hour,minute,second=xlrd.xldate_as_tuple(last_date_value,workbook.datemode)
        delta=datetime.date(year,month,day)-initial_date
        #put it in a dict
        data[delta.days]=worksheet.cell(row, value_column).value
        #print('world.indices['+str(delta.days)+']='+str(worksheet.cell(row, value_column).value)) #just for dengueme
    return data

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + datetime.timedelta(n)


def getValuesFromCsv(filename,start_date,end_date,value_column,verbose=True):
    with open(filename) as f:
        content = f.readlines()

    date_column=0
    date_format=''
    #guess the date format
    first_date=content[1].split(',')[date_column]
    if('-' in first_date):
        date_format="%Y-%m-%d"
    elif(':' in first_date):
        #date_format="%Y%m%d"
        print('Error data by the hour not supported!!!')
    else:
        date_format="%Y%m%d"

    #we build the value dict
    values={}
    content = [x.strip() for x in content[2:]]#strip trailing \n's and remove header
    for line in content:
        splitted_line=line.split(',')
        date=datetime.datetime.strptime(splitted_line[date_column], date_format).date()
        values[str(date)]=splitted_line[value_column]

    #we translate the dict to an ordered list(TODO:Find a better way!!!). Maybe collections.OrderedDict(sorted(values.items())
    values_list=[]
    for date in daterange(start_date,end_date):
        if(str(date) in values):
            value= values[str(date)]
            if(value==''):
                values_list.append(None)#no info
            else:
                values_list.append(float(value))
        else:
            values_list.append(None)#no info
            if(verbose): print ('No info for ' + str(date))
    return values_list

def getIndexesForPlot(filename,initial_date,date_column,value_column,filter_value=None,filter_column=None):
    indexes_data=getValuesFromXls(filename,initial_date,date_column,value_column,filter_value,filter_column)
    lists = sorted(indexes_data.items()) # sorted by key, return a list of tuples
    days, indexes = zip(*lists) # unpack a list of pairs into two tuples
    return np.asarray(days).tolist(),np.asarray(indexes).tolist()#to list

def getDatesTicksForPlot(start_date, time_range):
    locations = [d for d in time_range if (datetime.timedelta(days=d)+start_date).day==1]
    labels = [datetime.timedelta(days=d)+start_date for d in time_range if (datetime.timedelta(days=d)+start_date).day==1]
    return locations,labels

def  getMinTemperaturesFromCsv(filename,start_date,end_date):#in Kelvin
    return [x +273.15 if x is not None else None for x in getValuesFromCsv(filename,start_date,end_date,1)]

def  getAverageTemperaturesFromCsv(filename,start_date,end_date):#in Kelvin#TODO:change name to meanTemperature
    return [x +273.15 if x is not None else None for x in getValuesFromCsv(filename,start_date,end_date,2)]

def  getMaxTemperaturesFromCsv(filename,start_date,end_date):#in Kelvin
    return [x +273.15 if x is not None else None for x in getValuesFromCsv(filename,start_date,end_date,3)]

#convenience method that return an array with the min temperature at d+0.0 and the max temperature at d+0.5 and its time_domain
def  getMinMaxTemperaturesFromCsv(filename,start_date,end_date):#in Kelvin
    min_temperatures=getMinTemperaturesFromCsv(filename,start_date,end_date)
    max_temperatures=getMaxTemperaturesFromCsv(filename,start_date,end_date)
    days=(end_date - start_date).days
    time_domain=np.linspace(0,days-1,days*2 -1)
    return time_domain,[min_temperatures[int(d)] if d%1.==0. else max_temperatures[int(d)] for d in time_domain]

def  getPrecipitationsFromCsv(filename,start_date,end_date):#in mm
    return [x if not x<0 else None for x in getValuesFromCsv(filename,start_date,end_date,4)]

def  getRelativeHumidityFromCsv(filename,start_date,end_date):#in percentage
    return [x if x else None for x in getValuesFromCsv(filename,start_date,end_date,5)]

def  getMeanWindSpeedFromCsv(filename,start_date,end_date):#in km/h
    return [x if x else None for x in getValuesFromCsv(filename,start_date,end_date,7)]

def  getOvitrapEggsFromCsv(filename,start_date,end_date,ovitrap):#in km/h
    return [x for x in getValuesFromCsv(filename,start_date,end_date,ovitrap,False)]

def getDailyResults(time_range,RES,start_date,end_date):
    daily_RES=[]
    for d in range(0,(end_date-start_date).days):
        date_d=start_date+datetime.timedelta(days=d)
        mean_RES_d=RES[np.trunc(time_range)==d,:].mean(axis=0)#for the day "d", we calculate the daily mean of E,L,P,etc
        daily_RES.append(mean_RES_d)
    return np.array(daily_RES)

def saveResults(time_range,RES,start_date,end_date):
    filename='backup/previous_results/'+datetime.datetime.now().strftime('%Y-%m-%d__%H_%M_%S')+'.csv'
    file=open(filename,'w')
    daily_RES=getDailyResults(time_range,RES,start_date,end_date)
    for d,daily_RES_d in enumerate(daily_RES):
        date_d=start_date+datetime.timedelta(days=d)
        file.write(date_d.strftime('%Y-%m-%d')+','+','.join([str(value) for value in daily_RES_d ])+ '\n')
    return filename

def loadResults(filename,start_date):
    converters = {0: lambda d: (datetime.datetime.strptime(d, '%Y-%m-%d').date()-start_date).days}#date->days passed from start_date
    RES=np.loadtxt(filename,delimiter=',',converters=converters)
    return RES[:,1:]# 1: is to discard the date column

#
def getSurface(x=None,y=None,r=None):#surface in cm2. x,y,r must be in cm
    if(r):
        return (math.pi*r**2)
    elif(x and y):
        return x*y
    else:
        assert(False)

def getCapacity(x=None,y=None,r=None,z=None):#capacity in litres. x,y,z,r must be in cm
    if (r and z):
        return getSurface(r=r)*z * 1e-3 # cm3->litres
    elif(x and y and z):
        return getSurface(x=x,y=y)*z * 1e-3 # cm3->litres
    else:
        assert(False)
