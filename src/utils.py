from configparser import ConfigParser
from equations import f,diff_eqs
import matplotlib.dates
import collections
import numpy as np
import pylab as pl
import matplotlib
import datetime
import sys

###############################################################I/O###############################################################
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
    content = [x.strip() for x in content[1:]]#strip trailing \n's and remove header
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
    if(not os.path.isfile(filename)): return [0.0],[0]
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

def getLocations():
    config_parser = ConfigParser()
    config_parser.read('resources/get_weather.cfg')
    return config_parser.sections()

#Reduce the resolution to leave pixels as blocks of 100mx100m (10mx10m --> 100mx100m)
def getReducedMatrix(S):
    #Clip the image to leave rows and columns multiple of ten
    S=S[:S.shape[0]-S.shape[0]%10 , :S.shape[1]-S.shape[1]%10]
    n,m=int(S.shape[0]/10),int(S.shape[1]/10)
    B=np.array([np.hsplit(b,m) for b in  np.vsplit(S,n)])
    M=np.mean(B,axis=(2,3))
    return M

def getPreferenceMatrix():
    C=np.load('out/C.npy')

    #assign each class points. like S[C=2]=9,or S[C=5]=0 so class 2 is very good(grass or homes) we assign a ten. class 5 is very bad (cement)
    S=np.zeros(C.shape)
    #TODO: assign real scores!
    S[C==0]=1*0
    S[C==1]=2*0
    S[C==2]=3*0
    S[C==3]=4*1

    M=getReducedMatrix(S)#get a matrix of pixels 100mx100m (blocks)

    #Create the preference matrix
    P=np.zeros((M.shape[0],M.shape[1],8))
    P[:,:,0]=np.roll(M,1,axis=0)#up
    P[:,:,1]=np.roll(np.roll(M,-1,axis=1),1,axis=0)#up-right
    P[:,:,2]=np.roll(M,-1,axis=1)#right
    P[:,:,3]=np.roll(np.roll(M,-1,axis=1),-1,axis=0)# down-right
    P[:,:,4]=np.roll(M,-1,axis=0)#down
    P[:,:,5]=np.roll(np.roll(M,1,axis=1),-1,axis=0)#down-left
    P[:,:,6]=np.roll(M,1,axis=1)#left
    P[:,:,7]=np.roll(np.roll(M,1,axis=1),1,axis=0)#up-left
    PERPENDICULAR=(0,2,4,6)
    DIAGONAL=(1,3,5,7)
    P[:,:,PERPENDICULAR]=P[:,:,PERPENDICULAR]/np.maximum(1,np.sum(P[:,:,PERPENDICULAR],axis=2)[:,:,np.newaxis])*4.#normalize and multiply by 4
    P[:,:,DIAGONAL]=P[:,:,DIAGONAL]/np.maximum(1,np.sum(P[:,:,DIAGONAL],axis=2)[:,:,np.newaxis])*4.#normalize

    return P

###############################################################Equation Decorators###############################################################
class MetricsEquations:
    def __init__(self,model,diff_eqs):
        self.model=model
        self.diff_eqs=diff_eqs

    def __call__(self,Y,t,parameters):
        dY=self.diff_eqs(Y,t,parameters)
        time_range=self.model.time_range
        #account calls
        if(parameters.calls is None):
            parameters.calls=np.array([0]*len(time_range))
        parameters.calls[(np.abs(time_range-t)).argmin()]+=1
        #account negatives
        if(parameters.negatives is None):
            parameters.negatives=np.array([0]*len(time_range))
        if(np.any(Y<0)): parameters.negatives[(np.abs(time_range-t)).argmin()]+=1
        return dY

class ProgressEquations:
    def __init__(self,model,diff_eqs):
        self.model=model
        self.diff_eqs=diff_eqs
        self.t_max=np.max(model.time_range)
    def __call__(self,Y,t,parameters):
        sys.stdout.write("Completed: %d%%   \r" % ( t/self.t_max *100.) )
        sys.stdout.flush()
        return self.diff_eqs(Y,t,parameters)
###############################################################Plot################################################################
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
    safe_values=values.copy()
    safe_values[values==[None]]=np.nan#this is because in python3 None is unorderable, so values<0 breaks...
    values[safe_values<0]=0.#replace negatives with zeros, preserving Nones
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

def plot(model,subplots,plot_start_date,title='',figure=True,color=None):
    time_range=model.time_range
    RES=model.Y
    parameters=model.parameters
    T=parameters.weather.T
    p=parameters.weather.p
    RH=parameters.weather.RH
    vBS_mf,mf=parameters.vBS_mf,parameters.mf
    BS_a,vBS_h,vBS_s,vBS_d,n=parameters.BS_a,parameters.vBS_h,parameters.vBS_s,parameters.vBS_d,parameters.n
    EGG,LARVAE,PUPAE,ADULT1,ADULT2=parameters.EGG,parameters.LARVAE,parameters.PUPAE,parameters.ADULT1,parameters.ADULT2
    AEDIC_INDICES_FILENAME='data/private/Indices aedicos Historicos '+parameters.location['name']+'.xlsx'

    if(figure): pl.figure()
    pl.subplots_adjust(top=0.95,hspace=0.28)
    ax1=None
    for i,subplot in enumerate(subplots):
        subplot_id=len(subplots)*100 + 10 + (i+1)
        if(i==0):
            ax1=pl.subplot(subplot_id)
            date_range=[datetime.timedelta(days=d)+datetime.datetime.combine(model.start_date,datetime.time()) for d in time_range]
            ax1.xaxis.set_major_locator( matplotlib.dates.MonthLocator() )
            ax1.xaxis.set_major_formatter( matplotlib.dates.DateFormatter('%Y-%m-%d') )
            ax1.xaxis.set_minor_locator( matplotlib.dates.WeekdayLocator(byweekday=matplotlib.dates.MO) )
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
        if ('A1+A2' in subplot): pl.plot(date_range,applyFs(RES[:,ADULT2]+RES[:,ADULT1],subplot), '-m',color=color, label='A1+A2'+title)

        #derivate
        dY=np.zeros(RES.shape)
        if(np.any([a==b for a in subplot for b in ['dE','dL','dP','dA1','dA1','dW']])):#if any of dE,...,dA2 is asked, calculate the whole dY
            for i,t in enumerate(time_range):
                dY[i,:]=diff_eqs(RES[i,:],t,model.parameters)

        if ('dE' in subplot): pl.plot(date_range,applyFs(dY[:,EGG],subplot), '-k', label='dE')
        if ('dL' in subplot): pl.plot(date_range,applyFs(dY[:,LARVAE],subplot), '-r', label='dL')
        if ('dP' in subplot): pl.plot(date_range,applyFs(dY[:,PUPAE],subplot), '-g', label='dP')
        if ('dA1' in subplot): pl.plot(date_range,applyFs(dY[:,ADULT1],subplot), '-b', label='dA1')
        if ('dA2' in subplot): pl.plot(date_range,applyFs(dY[:,ADULT2],subplot), '-m', label='dA2')
        if ('dW' in subplot): pl.plot(date_range,applyFs(dY[:,WATER],subplot), '-b', label='dW')

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
        if('lwL' in subplot):
            lwL=np.array([RES[(np.abs(time_range-t)).argmin(),LARVAE]-RES[(np.abs(time_range-(t-7))).argmin(),LARVAE] for t in time_range])
            pl.plot(date_range, applyFs(lwL,subplot), '-', label='L(t)-L(t-7)')
        pl.ylabel('')

        #Complete lifecycle
        if('clc' in subplot):
            pl.plot(date_range,[sum([1./model.R_D(stage,model.T(t)) for stage in [EGG,LARVAE,PUPAE,ADULT1,ADULT2]]) for  t in time_range],label='Complete life cicle(from being an egg to the second oviposition)')
            pl.ylabel('days')
        #Water in containers(in L)
        if ('W' in subplot):
            mW=parameters.vW(time_range)
            pl.plot(date_range,applyFs(mW,subplot), label='W(t)')
            pl.ylabel('cm.')

        #manually_filled(in mm.)
        if ('mf' in subplot):
            pl.plot(date_range,applyFs(np.array([mf(t)*vBS_mf*vBS_h*10. for t in time_range]),subplot), label='mf(t)')
            pl.ylabel('mm./day')

        #spaa vs cimsim
        if ('spaavscimsim' in subplot):
            for i in range(0,n):
                pl.plot(time_range,RES[:,WATER[i] ], label='W(t) for %scm, %scm^2, %s%%'%(vBS_h[i],vBS_s[i],vBS_d[i]*100.) )#L->ml->mm->cm
            pl.plot(getValuesFromCsv('data/test/cimsim_containers_2015_se09.csv',model.start_date,model.end_date,1,verbose=False),label='CIMSiM')

        #Temperature in K
        if ('T' in subplot):
            pl.plot(date_range,applyFs(np.array([T(t) for t in time_range]),subplot),color=color, label='Temperature'+title)
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
            pl.plot(date_range,[f( RES[(np.abs(time_range-t)).argmin(),WATER] , vBS_d ) for t in time_range], label='f(vW,vBS_d)')
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
        pl.title(title)

def showPlot():
    return pl.show()

###############################################################Animation################################################################
from PIL import ImageFont, ImageDraw,Image
from moviepy.video.io.bindings import PIL_to_npimage
import moviepy.editor as mpy
#from https://gist.github.com/Zulko/06f49f075fd00e99b4e6#file-moviepy_time_accuracy-py-L33-L39
def addText(matrix,text):
    im = Image.fromarray( np.uint8(matrix))
    draw = ImageDraw.Draw(im)
    draw.text((50, 25), str(text))
    return PIL_to_npimage(im)

def createAnimation(out_filename,matrix,getTitle,duration):
    R=np.load('out/R.npy')#Load the original Raster
    R=np.moveaxis(R,0,-1)#(4,n,m) ----> (n,m,4)
    R=getReducedMatrix(R[:,:,0:3])#10mx10m ----> 100mx100m
    R=np.clip(R/(2*R.mean()),0,1)#this is just to make the base image look nice.
    red = np.array([1,0,0]).transpose()
    def makeFrame(t):
        frame=255*(0.5*R + 0.5*red*matrix[int(t),:,:,np.newaxis])
        return addText(frame, getTitle(int(t)))

    animation = mpy.VideoClip(makeFrame, duration=duration)
    animation.write_videofile(out_filename+'.mp4', fps=15)
