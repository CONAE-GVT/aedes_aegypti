from configparser import ConfigParser
import plotly.graph_objs as go
from equations import diff_eqs,ovsp,vR_D
import plotly.offline as ply
import collections
import numpy as np
import datetime
import tempfile
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

def getOvitrapEggsFromCsv2(filename,start_date,end_date,column):#amount
    lines=[line.strip().split(',') for line in open(filename).readlines()[1:]]
    values={}
    last_date=None
    for line in lines:
        if(line[0]):
            last_date=date=datetime.datetime.strptime(line[0], '%Y-%m-%d').date()
        else:
            date=last_date

        if(line[column]):
            value=float(line[column])
        else:
            value=None

        if(date in values):
            values[date].append(value)
        else:
            values[date]=[value]
    return values

def getStartEndDates(filename):
    dates=[line.split(',')[0] for line in open(filename,'r').readlines()]
    return datetime.datetime.strptime(dates[1], '%Y-%m-%d').date(),datetime.datetime.strptime(dates[-1], '%Y-%m-%d').date()

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
def getReducedMatrix(S,block_size=10):
    #Clip the image to leave rows and columns multiple of ten
    S=S[:S.shape[0]-S.shape[0]%block_size , :S.shape[1]-S.shape[1]%block_size]
    n,m=int(S.shape[0]/block_size),int(S.shape[1]/block_size)
    B=np.array([np.hsplit(b,m) for b in  np.vsplit(S,n)])
    M=np.mean(B,axis=(2,3))
    return M

import gdal
def getPreferenceMatrix(raster_filename,patch_size):
    raster=gdal.Open(raster_filename, gdal.GA_ReadOnly)
    pixel_size_X,pixel_size_Y= raster.GetGeoTransform()[1], -raster.GetGeoTransform()[5]
    assert pixel_size_X==pixel_size_Y#right now we just tested on images with squared pixels, but it shouldn't be hard to accepts different values
    S=raster.ReadAsArray()
    S[S<-1e30]=0.#no data --> 0
    block_size=int(round(patch_size/int(pixel_size_X)))
    warning=None
    if(block_size*pixel_size_X != patch_size): warning='Using a patch of %smx%sm instead of %smx%sm'%(int(block_size*pixel_size_X),int(block_size*pixel_size_X),patch_size,patch_size)
    M=getReducedMatrix(S,block_size=int(round(patch_size/int(pixel_size_X))))#get a matrix of pixels 100mx100m (blocks)

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

    return P,warning

from skimage import transform
def getY0FactorMatrix(height,width):
    Y0_f=np.load('out/Y0_f.npy')
    return transform.resize(Y0_f,(height,width),mode='constant')

###############################################################Equation Decorators###############################################################
class MetricsEquations:
    def __init__(self,model,diff_eqs):
        self.model=model
        self.diff_eqs=diff_eqs
        model.parameters.calls=None
        model.parameters.negatives=None

    def __call__(self,Y,t,h,parameters):
        dY=self.diff_eqs(Y,t,h,parameters)
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

class EgnCorrector:
    def __init__(self,p,BS_a,start_date,end_date):
        self.ovitrap_eggs=np.array(getOvitrapEggsFromCsv('data/private/ovitrampas_2017-2018.csv',start_date,end_date,153))
        self.p=p
        self.BS_a=BS_a
        self.egn_corrections=[]

    def __call__(self,egn,cycle_A,t):#TODO:take into account amount of containers
        if(self.p==0): return egn
        OE_week=self.ovitrap_eggs[int(t):int(t)+7]#TODO:check if int(t) is ok. maybe it's ceil.#also use an spline and derivative instead of /7
        OE_week=OE_week[OE_week!=[None]]
        if(len(OE_week)==0): return 19.#TODO:what should we return?
        p,BS_A=self.p,self.BS_a
        dOEdt=np.asscalar(OE_week)/7. * self.BS_a
        self.egn_corrections.append(dOEdt/cycle_A)
        p=self.p
        return (dOEdt/cycle_A)*p + (1-p)*egn

class OEquations:
    def __init__(self,model,diff_eqs):
        self.model=model
        self.diff_eqs=diff_eqs
        m,n=model.parameters.m,model.parameters.n
        OVIPOSITION_START=len(model.parameters.initial_condition)
        model.parameters.OVIPOSITION=slice(OVIPOSITION_START,OVIPOSITION_START + m*n)#in R^(mxn)
        model.parameters.initial_condition=np.append(model.parameters.initial_condition,np.zeros(m*n))
        assert model.parameters.initial_condition[model.parameters.ADULT1]==model.parameters.initial_condition[model.parameters.ADULT2]==0#otherwise the above statement is false.

    def __call__(self,Y,t,h,parameters):
        dY=self.diff_eqs(Y,t,h,parameters)
        T_t=parameters.weather.T(t)
        elr,lpr,par,ovr1,ovr2=vR_D(T_t)
        BS_lh,vBS_d,m,n,mBS_l=parameters.BS_lh,parameters.vBS_d,parameters.m,parameters.n,parameters.mBS_l
        ADULT1,ADULT2,WATER,OVIPOSITION=parameters.ADULT1,parameters.ADULT2,parameters.WATER,parameters.OVIPOSITION

        A1,A2,vW_t=Y[ADULT1],Y[ADULT2],Y[WATER]
        vW_l=vW_t/BS_lh

        egn=63.0
        ovsp_t=ovsp(vW_t,vBS_d,vW_l,mBS_l)
        egn_c=parameters.egnCorrector(egn,ovr1 *A1  + ovr2* A2, t)
        dO=egn_c*( ovr1 *A1  + ovr2* A2)*ovsp_t
        dY[OVIPOSITION]=dO.transpose().reshape((1,m*n))

        return dY

###############################################################Plot################################################################
def normalize(values):#TODO:change name.
    return (values-values.min())/(values.max()-values.min())

def safeNormalize(values):#take care of None values
    values=np.array([v/(values[values!=[None]].max()) if v!=None else None for v in values])
    return values

def safeAdd(values):
    if(len(values.shape)==3): return np.sum(values,axis=(1,2))
    if(len(values.shape)!=2):
        return values
    else:
        return np.sum(values,axis=1)

def replaceNegativesWithZeros(values):
    safe_values=values.copy()
    safe_values[values==[None]]=np.nan#this is because in python3 None is unorderable, so values<0 breaks...
    values[safe_values<0]=0.#replace negatives with zeros, preserving Nones
    return values

def noneMean(values):
    values=np.array(values)
    values=values[values!=[None]]
    if(len(values)!=0):
        return np.mean(values)
    else:
        return None

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

def plot(model,subplots,plot_start_date=None,color=None):
    time_range=model.time_range
    RES=model.Y
    parameters=model.parameters
    T=parameters.weather.T
    p=parameters.weather.p
    RH=parameters.weather.RH
    vBS_mf,mf=parameters.vBS_mf,parameters.mf
    BS_a,BS_lh,vBS_h,vBS_s,vBS_d,m,n=parameters.BS_a,parameters.BS_lh,parameters.vBS_h,parameters.vBS_s,parameters.vBS_d,parameters.m,parameters.n
    EGG,LARVAE,PUPAE,ADULT1,ADULT2,WATER=parameters.EGG,parameters.LARVAE,parameters.PUPAE,parameters.ADULT1,parameters.ADULT2,parameters.WATER
    data=[]

    date_range=[datetime.timedelta(days=d)+datetime.datetime.combine(model.start_date,datetime.time()) for d in time_range]
    for i,subplot in enumerate(subplots):
        if(plot_start_date):
            time_range,RES,date_range=subData(time_range,RES,date_range,plot_start_date)

        #Amount of larvaes,pupaes and adults
        if ('E_i' in subplot):
            for i,y in enumerate(applyFs(RES[:,EGG],subplot).transpose()):
                bs_i=i%m
                label='E in [%.1f,%.1f)cm'%(bs_i*BS_lh,(bs_i+1)*BS_lh)
                data.append(go.Bar(x=date_range,y=y,name=label,text=label))
                #data.append(go.Scatter(x=date_range,y=y, name='E in [%.1f,%.1f)cm'%(bs_i*BS_lh,(bs_i+1)*BS_lh)))
        if ('E' in subplot): data.append(go.Scatter(x=date_range,y=applyFs(RES[:,EGG],subplot)*1e-2, name='Eggs . 10e-2'))
        if ('L' in subplot): data.append(go.Scatter(x=date_range,y=applyFs(RES[:,LARVAE],subplot), name='L'))
        if ('P' in subplot): data.append(go.Scatter(x=date_range,y=applyFs(RES[:,PUPAE],subplot), name='P'))
        if ('A1' in subplot): data.append(go.Scatter(x=date_range,y=applyFs(RES[:,ADULT1],subplot), name='A1'))
        if ('A2' in subplot): data.append(go.Scatter(x=date_range,y=applyFs(RES[:,ADULT2],subplot), name='A2'))
        if ('A1+A2' in subplot): data.append(go.Scatter(x=date_range,y=applyFs(RES[:,ADULT2]+RES[:,ADULT1],subplot), name=subplot['A1+A2'] or 'Adults',line = dict(color= (color)) ))

        if('lwO' in subplot):
            Y=RES#it should be Y everywhere....
            indexOf=lambda t: (np.abs(time_range-t)).argmin()
            OVIPOSITION=model.parameters.OVIPOSITION
            O=Y[:,OVIPOSITION]
            lwO=np.array([Y[indexOf(t),OVIPOSITION]-Y[indexOf(t-7),OVIPOSITION] for t in time_range])/BS_a
            lwO_mean=np.array([lwO[indexOf(t-7):indexOf(t+7)].mean(axis=0) for t in time_range])
            lwO_std =np.array([lwO[indexOf(t-7):indexOf(t+7)].std(axis=0) for t in time_range])
            data.append(go.Scatter(x=date_range, y=applyFs(lwO,subplot), name=subplot['lwO'] or 'lwO',line = dict(color= (color))  ))

        #delta Eggs
        if('lwE' in subplot):
            lwE=np.array([RES[(np.abs(time_range-t)).argmin(),EGG]-RES[(np.abs(time_range-(t-7))).argmin(),EGG] for t in time_range])
            lwE_mean=np.array([lwE[(np.abs(time_range-(t-7))).argmin():(np.abs(time_range-(t+7))).argmin()].mean(axis=0) for t in time_range])/BS_a
            lwE_std =np.array([lwE[(np.abs(time_range-(t-7))).argmin():(np.abs(time_range-(t+7))).argmin()].std(axis=0) for t in time_range])/BS_a
            data.append(go.Scatter(x=date_range, y=applyFs(lwE_mean,subplot), name='E(t)-E(t-7) mean'))

        #Complete lifecycle
        if('clc' in subplot):
            data.append(go.Scatter(x=date_range,y=[sum([1./model.R_D(stage,model.T(t)) for stage in [EGG,LARVAE,PUPAE,ADULT1,ADULT2]]) for  t in time_range],label='Complete life cicle(from being an egg to the second oviposition) in days'))
        #Water in containers(in L)
        if ('W' in subplot):
            mW=RES[:,WATER]
            for y in applyFs(mW,subplot).transpose():
                data.append(go.Scatter(x=date_range,y=y, name=subplot['W'] or 'W(t) in cm.',line = dict(color= (color)) ))

        #Temperature in K
        if ('T' in subplot):
            data.append(go.Scatter(x=date_range,y=applyFs(np.array([T(t)- 273.15 for t in time_range]),subplot),name='Temperature in C'))

        #precipitations rate(in mm./day)
        if ('p' in subplot):
            data.append(go.Scatter(x=date_range,y=applyFs(np.array([p(t+0.5) for t in time_range]),subplot), name='p(t+1/2) in mm./day'))


        #precipitations accumulated(in mm.)
        if ('pa' in subplot):
            precipitations = getPrecipitationsFromCsv('data/public/'+model.parameters.location['name']+'.csv',date_range[0].date(),date_range[-1].date())
            data.append(go.Bar(x=date_range,y=precipitations,name='Accumulated precipitations in mm.'))

        #relative Humidity
        if ('RH' in subplot):
            data.append(go.Scatter(x=date_range,y=applyFs(np.array([RH(t) for t in time_range]),subplot), label='RH(t) in %'))

        if('O' in subplot):
            for ovitrap_id in subplot['O']:
                values=getOvitrapEggsFromCsv2('data/private/ovitrampas_2017-2018.full.csv',model.start_date,model.end_date,ovitrap_id)
                ovitrap_dates=np.array([k for k in values.keys()])
                ovi=np.array([noneMean(values[date]) for date in ovitrap_dates])
                data.append(go.Scatter(x=ovitrap_dates[ovi!=[None]], y=applyFs(ovi,subplot)[ovi!=[None]], name='Ovitrap %s eggs'%ovitrap_id, mode = 'markers'))

        if('cd' in subplot):#current date
            location=model.parameters.location['name']
            if('cordoba.full.weather-' in location):
                current_date=datetime.datetime.strptime(location.replace('cordoba.full.weather-',''),'%Y-%m-%d').date()
            else:
                current_date=datetime.date.today()
            data.append(go.Scatter(x=[current_date,current_date],y=[0,600],name='Current Date',mode='lines'));
        #debugging plots
        #Calls
        if ('c' in subplot):
            data.append(go.Scatter(xdate_range,y=model.parameters.calls, label='calls'))
            print('Calls: %s'%sum(model.parameters.calls))
        #Negatives
        if ('n' in subplot):
            data.append(go.Scatter(x=date_range,y=model.parameters.negatives, label='negatives'))
            print('Negatives: %s'%sum(model.parameters.negatives))

    return  data

def showPlot(data,title='',xaxis_title='',yaxis_title=''):
    layout=go.Layout(title=title,barmode='stack',
                    xaxis = dict(title = xaxis_title),
                    yaxis = dict(title = yaxis_title))
    ply.plot(go.Figure(data=data,layout=layout), filename=tempfile.NamedTemporaryFile(prefix='plot_').name)

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
    S=gdal.Open('data/public/goodness/predict.backgr.out_mean.tif', gdal.GA_ReadOnly).ReadAsArray()
    S[S<-1e30]=0.#no data --> 0
    R=getReducedMatrix(S,block_size=int(round(100/6)))#get a matrix of pixels 100mx100m (blocks)
    red = np.array([1,0,0]).transpose()
    grey= (np.array([119,136,153])/255).transpose()
    def makeFrame(t):
        frame=255*(0.5*grey*R[:,:,np.newaxis] + 0.5*red*matrix[int(t),:,:,np.newaxis])
        return addText(frame, getTitle(int(t)))

    animation = mpy.VideoClip(makeFrame, duration=duration)
    animation.write_videofile(out_filename+'.mp4', fps=15)



###############################################################Metrics#####################################################################
def rmse(x,y):
    assert len(x)==len(y)
    return delta_e(x,y)/np.sqrt(len(x))#Barnston, A., (1992). “Correspondence among the Correlation [root mean square error] and Heidke Verification Measures; Refinement of the Heidke Score.” Notes and Correspondence, Climate Analysis Center.

#Astrain-Ojeda 2017 "Medidas de disimilitud en series temporales"
#Chouakria-Nagabhushan 2007 "Adaptive dissimilarity index for measuring time series proximity"
def cort(x,y):
    assert len(x)==len(y)
    return np.sum((x[:-1]-x[1:])*(y[:-1]-y[1:]))/(  np.sqrt(np.sum((x[:-1]-x[1:])**2)) * np.sqrt(np.sum((y[:-1]-y[1:])**2)) )

def delta_e(x,y):
    assert len(x)==len(y)
    return np.sqrt( np.sum((x-y)**2) )

def D(x,y):
    return f(cort(x,y))*delta_e(x,y)

def f(x,k=3):
    return 2/(1+np.exp(k*x))
