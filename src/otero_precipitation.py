#coding: utf-8
import scipy.integrate as spi
from scipy import interpolate
from config import Configuration
from weather import Weather
from bunch import Bunch
import numpy as np
import spatial_equations as equations
import datetime
import utils
import rk

class Model:
    def __init__(self, configuration=Configuration('resources/otero_precipitation.cfg')):
        self.configuration=configuration

        self.parameters=Bunch()
        self.parameters.BS_a=configuration.getFloat('breeding_site','amount')
        self.parameters.vBS_h=configuration.getArray('breeding_site','height')#in cm
        self.parameters.vBS_s=configuration.getArray('breeding_site','surface')#in cm^2
        self.parameters.vBS_d=configuration.getArray('breeding_site','distribution')#distribution of BS. Sum must be equals to 1
        self.parameters.vBS_W0=configuration.getArray('breeding_site','initial_water')
        self.parameters.vBS_mf=configuration.getArray('breeding_site','manually_filled')#in percentage of capacity
        self.parameters.n=len(self.parameters.vBS_d)

        n=self.parameters.n
        self.parameters.EGG=slice(0,n)#in R^n
        self.parameters.LARVAE=slice(n,2*n)#in R^n
        self.parameters.PUPAE=slice(2*n,3*n)#in R^n
        self.parameters.ADULT1=3*n#in R
        self.parameters.FLYER=3*n+1#in R
        self.parameters.ADULT2=3*n+2#in R
        self.parameters.vAlpha0=configuration.getArray('biology','alpha0')#constant to be fitted

        #Cordoba
        self.parameters.location={'name':configuration.getString('location','name')}
        self.start_date=configuration.getDate('simulation','start_date')
        self.end_date=configuration.getDate('simulation','end_date')
        self.time_range = np.linspace(0, (self.end_date - self.start_date).days-1, (self.end_date - self.start_date).days )
        self.parameters.initial_condition=configuration.getArray('simulation','initial_condition')

        WEATHER_DATA_FILENAME='data/public/'+self.parameters.location['name']+'.csv'
        self.parameters.weather=Weather(WEATHER_DATA_FILENAME,self.start_date,self.end_date)

        self.parameters.mf=self.parameters.weather.getAsLambdaFunction(self.parameters.weather.aps, [0,0,0,0,0,0,1.]* int( (self.end_date - self.start_date).days/7 +1) )
        W = spi.odeint(equations.waterEquations,self.parameters.vBS_W0,self.time_range,hmax=1.0,args=(self.parameters,))
        self.parameters.vW=interpolate.interp1d(self.time_range,W,axis=0)

        self.parameters.P=utils.getPreferenceMatrix()
        self.validate()

    def validate(self):
        mean_temperatures=np.array([self.parameters.weather.T(t) for t in self.time_range])
        lower_bound=mean_temperatures[mean_temperatures<278.]
        upper_bound=mean_temperatures[mean_temperatures>303.]
        if(lower_bound.size>0 or upper_bound.size>0):
            print('# WARNING: Temperature out of model\'s valid range:T<278:%s T>303:%s'%(lower_bound.size,upper_bound.size))

    def save(self):
        #save results
        results_filename='data/test/previous_results/'+self.configuration.getString('location','name')+'-'+datetime.datetime.now().strftime('%Y-%m-%d__%H_%M_%S')+'.csv'
        file=open(results_filename,'w')
        daily_Y=utils.getDailyResults(self.time_range,self.Y,self.start_date,self.end_date)
        for d,daily_Y_d in enumerate(daily_Y):
            date_d=self.start_date+datetime.timedelta(days=d)
            file.write(date_d.strftime('%Y-%m-%d')+','+','.join([str(value) for value in daily_Y_d ])+ '\n')
        #save config
        self.configuration.save(results_filename.replace('.csv','.cfg'))

        return results_filename


    def solveEquations(self,equations=equations.diff_eqs,method='odeint'):
        time_range=self.time_range
        n=self.parameters.n

        HEIGHT,WIDTH=self.parameters.P.shape[:2]
        tmp=np.zeros((HEIGHT,WIDTH,3*n + 3))
        tmp[int(HEIGHT/2),int(WIDTH/2),:]=1.
        initial_condition=(self.parameters.initial_condition*tmp).reshape((HEIGHT*WIDTH*(3*n + 3) ))#TODO:check that this does what we expect.
        self.parameters.vBS_a=self.parameters.BS_a*np.ones((HEIGHT,WIDTH))#np.random.random((WIDTH,HEIGHT))#TODO:do something about this...
        self.parameters.vBS_d=self.parameters.vBS_d*np.ones((HEIGHT,WIDTH,n))
        self.parameters.vAlpha=(self.parameters.vAlpha0*np.ones((HEIGHT,WIDTH,n)) )/self.parameters.vBS_a[:,:,np.newaxis]
        self.parameters.ovr=np.ones((HEIGHT,WIDTH))/0.229#TODO:implement!!!!
        Y=None

        if(method=='odeint'):
            Y = spi.odeint(equations,initial_condition,time_range,hmax=1.0,args=(self.parameters,))#the ',' in (parameters,) is very important! '(parameters)' or tuple(parameters) doesn't work#TODO: this is because it calls aps out of it's domain.Find a better way.
        elif(method=='rk'):
            Y=rk.solve(equations,initial_condition,time_range,args=(self.parameters,),steps=20)
        elif(method=='dopri'):
            Y=rk.scipy_solve(equations,initial_condition,time_range,'dopri',{'max_step':time_range[1]-time_range[0],'rtol':1e-3, 'atol':1e-6}, args=(self.parameters,))

        self.Y=Y
        return time_range,initial_condition,Y
