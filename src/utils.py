from equations import f,diff_eqs
import matplotlib.dates
import collections
import numpy as np
import pylab as pl
import matplotlib
import datetime
import math
import xlrd


###############################################################I/O###############################################################
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

def  getOvitrapEggsFromCsv(filename,start_date,end_date,ovitrap):#amount
    return [x for x in getValuesFromCsv(filename,start_date,end_date,ovitrap,False)]

def getDailyResults(time_range,RES,start_date,end_date):
    daily_RES=[]
    for d in range(0,(end_date-start_date).days):
        date_d=start_date+datetime.timedelta(days=d)
        mean_RES_d=RES[np.trunc(time_range)==d,:].mean(axis=0)#for the day "d", we calculate the daily mean of E,L,P,etc
        daily_RES.append(mean_RES_d)
    return np.array(daily_RES)

def loadResults(filename,start_date):
    converters = {0: lambda d: (datetime.datetime.strptime(d.decode('ascii'), '%Y-%m-%d').date()-start_date).days}#date->days passed from start_date
    RES=np.loadtxt(filename,delimiter=',',converters=converters)
    return RES[:,1:]# 1: is to discard the date column

###############################################################Misc.###############################################################
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

###############################################################Plot###############################################################
def normalize(values):#TODO:change name.
    return (values-values.min())/(values.max()-values.min())

def safeNormalize(values):#take care of None values
    values=np.array([v/(values[values!=[None]].max()) if v!=None else None for v in values])
    return values

def safeAdd(values):
    if(len(values.shape)!=2):
        return values
    else:
        return np.sum(values,axis=1)

def replaceNegativesWithZeros(values):
    values[np.logical_and(values<0,values!=[None])]=0.#replace negatives with zeros, preserving Nones
    return values

def applyFs(values,subplot):
    if(type(subplot) is dict): subplot=subplot.values()#if we pass a dict, we need to iterate over the values(were the functions are.)
    for f_list in subplot:#first we assume all elements
        for f in f_list:   #in subplots are a list of f to apply
            if(callable(f)):# and then we check if it's actually a function
                values=f(values)
    return values

def subData(time_range,Y,date_range,an_start_date):
    #get the index of an_start_date
    index=None
    for i,date in enumerate(date_range):
        if(date.date()==an_start_date):
            index=i
            break#conserve the first one.
    return time_range[index:],Y[index:,:],date_range[index:]

def plot(model,subplots,plot_start_date):
    time_range=model.time_range
    RES=model.Y
    parameters=model.parameters
    T=parameters.weather.T
    p=parameters.weather.p
    RH=parameters.weather.RH
    BS_a,vBS_oc,vBS_ic,vBS_od,vBS_id,vBS_os,n,m=parameters.BS_a,parameters.vBS_oc,parameters.vBS_ic,parameters.vBS_od,parameters.vBS_id,parameters.vBS_os,parameters.n,parameters.m
    EGG,LARVAE,PUPAE,ADULT1,ADULT2,WATER=parameters.EGG,parameters.LARVAE,parameters.PUPAE,parameters.ADULT1,parameters.ADULT2,parameters.WATER
    AEDIC_INDICES_FILENAME='data/private/Indices aedicos Historicos '+parameters.location['name']+'.xlsx'

    pl.figure()
    pl.subplots_adjust(top=0.95,hspace=0.28)
    ax1=None
    for i,subplot in enumerate(subplots):
        subplot_id=len(subplots)*100 + 10 + (i+1)
        if(i==0):
            ax1=pl.subplot(subplot_id)
            date_range=[datetime.timedelta(days=d)+datetime.datetime.combine(model.start_date,datetime.time()) for d in time_range]
            ax1.xaxis.set_major_locator( matplotlib.dates.MonthLocator() )
            ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%Y-%m-%d') )
            ax1.xaxis.set_minor_locator( matplotlib.dates.DayLocator() )
        else:
            pl.subplot(subplot_id,sharex=ax1)#sharex to zoom all subplots if one is zoomed

        if(plot_start_date):
            time_range,RES,date_range=subData(time_range,RES,date_range,plot_start_date)

        #Amount of larvaes,pupaes and adults
        if ('E' in subplot): pl.plot(date_range,applyFs(RES[:,EGG],subplot), '-k', label='E')
        if ('L' in subplot): pl.plot(date_range,applyFs(RES[:,LARVAE],subplot), '-r', label='L')
        if ('P' in subplot): pl.plot(date_range,applyFs(RES[:,PUPAE],subplot), '-g', label='P')
        if ('A1' in subplot): pl.plot(date_range,applyFs(RES[:,ADULT1],subplot), '-b', label='A1')
        if ('A2' in subplot): pl.plot(date_range,applyFs(RES[:,ADULT2],subplot), '-m', label='A2')
        if ('A1+A2' in subplot): pl.plot(date_range,applyFs(RES[:,ADULT2]+RES[:,ADULT1],subplot), '-m', label='A1+A2')

        #derivate
        dY=np.zeros(RES.shape)
        if(np.any([a==b for a in subplot for b in ['dE','dL','dP','dA1','dA1']])):#if any of dE,...,dA2 is asked, calculate the whole dY
            for i,t in enumerate(time_range):
                dY[i,:]=diff_eqs(RES[i,:],t,model.parameters)

        if ('dE' in subplot): pl.plot(date_range,applyFs(dY[:,EGG],subplot), '-k', label='dE')
        if ('dL' in subplot): pl.plot(date_range,applyFs(dY[:,LARVAE],subplot), '-r', label='dL')
        if ('dP' in subplot): pl.plot(date_range,applyFs(dY[:,PUPAE],subplot), '-g', label='dP')
        if ('dA1' in subplot): pl.plot(date_range,applyFs(dY[:,ADULT1],subplot), '-b', label='dA1')
        if ('dA2' in subplot): pl.plot(date_range,applyFs(dY[:,ADULT2],subplot), '-m', label='dA2')

        #larvae indices
        if ('LI' in subplot):
            ri_days, ris=np.array(getIndexesForPlot(AEDIC_INDICES_FILENAME,model.start_date,0,5))
            pl.plot([datetime.timedelta(days=d)+datetime.datetime.combine(model.start_date,datetime.time()) for d in ri_days], applyFs(ris,subplot), '^y', label='Recip. Indices',clip_on=False, zorder=100,markersize=8)

        #ovitraps
        if('O' in subplot):
            for i in subplot['O']:
                ovitrap_eggs=np.array(getOvitrapEggsFromCsv('data/private/Datos sensores de oviposicion.NO.csv',model.start_date,model.end_date,i))
                pl.plot([datetime.timedelta(days=d)+datetime.datetime.combine(model.start_date,datetime.time()) for d in range(0,len(ovitrap_eggs))], applyFs(ovitrap_eggs,subplot), '^', label='Ovitrap %s eggs'%i,clip_on=False, zorder=100,markersize=8)

        #delta Eggs
        if('lwE' in subplot):
            lwE=np.array([RES[(np.abs(time_range-t)).argmin(),EGG]-RES[(np.abs(time_range-(t-7))).argmin(),EGG] for t in time_range])
            pl.plot(date_range, applyFs(lwE,subplot), '-', label='E(t)-E(t-7)')
        pl.ylabel('')

        #Complete lifecycle
        if('clc' in subplot):
            pl.plot(date_range,[sum([1./model.R_D(stage,model.T(t)) for stage in [EGG,LARVAE,PUPAE,ADULT1,ADULT2]]) for  t in time_range],label='Complete life cicle(from being an egg to the second oviposition)')
            pl.ylabel('days')
        #Water in containers(in L)
        if ('W' in subplot):
            for i in range(0,n):
                pl.plot(date_range,applyFs(RES[:,WATER[i] ],subplot), label='W(t) for %sL, %scm^2, %s%%'%(vBS_oc[i],vBS_os[i],vBS_od[i]*100.) )
            pl.ylabel('Litres')

        #spaa vs cimsim
        if ('spaavscimsim' in subplot):
            for i in range(0,n):
                pl.plot(time_range,RES[:,WATER[i] ]*1000.0/vBS_os[i], label='W(t) for %sL, %scm^2, %s%%'%(vBS_oc[i],vBS_os[i],vBS_od[i]*100.) )#L->ml->mm->cm
            pl.plot(getValuesFromCsv('data/test/cimsim_containers_2015_se09.csv',model.start_date,model.end_date,1,verbose=False),label='CIMSiM')

        #Temperature in K
        if ('T' in subplot):
            pl.plot(date_range,applyFs(np.array([T(t) for t in time_range]),subplot), label='Temperature')
            pl.ylabel('K')

        #precipitations(in mm.)
        if ('p' in subplot):
            pl.plot(date_range,applyFs(np.array([p(t) for t in time_range]),subplot),'-b', label='p(t)')
            pl.ylabel('mm./day')

        #relative Humidity
        if ('RH' in subplot):
            pl.plot(date_range,applyFs(np.array([RH(t) for t in time_range]),subplot), label='RH(t)')
            pl.ylabel('%')

        #f
        if ('b' in subplot):
            pl.plot(date_range,[f( np.concatenate((RES[(np.abs(time_range-t)).argmin(),WATER],vBS_ic)), np.concatenate((vBS_od,vBS_id)) ) for t in time_range], label='f(vW,vBS_d)')
            pl.ylabel('')

        #debugging plots
        #Calls
        if ('c' in subplot):
            pl.plot(date_range,model.parameters.calls, label='calls')
            pl.ylabel('')
            print('Calls: %s'%sum(model.parameters.calls))
        #Negatives
        if ('n' in subplot):
            pl.plot(date_range,model.parameters.negatives, label='negatives')
            pl.ylabel('')
            print('Negatives: %s'%sum(model.parameters.negatives))

        #common to all subplots
        pl.xlabel('Time(in days starting in July)')
        pl.legend(loc=0)
        pl.xticks(rotation='vertical')

from mpl_toolkits import mplot3d
def plot3D(xline,yline,zline):
    ax = pl.axes(projection='3d')
    pl.xlabel('time')
    pl.ylabel('W')
    #pl.zlabel('E')
    ax.plot3D(xline, yline, zline, 'gray')


def showPlot():
    return pl.show()
