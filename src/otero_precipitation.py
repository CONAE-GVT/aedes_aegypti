#coding: utf-8
import scipy.integrate as spi
from config import Configuration
from weather import Weather
from bunch import Bunch
import numpy as np
import equations
import datetime
import utils
import rk

class Model:
    def __init__(self, configuration=Configuration('resources/otero_precipitation.cfg')):
        self.parameters=Bunch()
        self.parameters.BS_a=configuration.getFloat('breeding_site','amount')
        self.parameters.vBS_oc=configuration.getArray('breeding_site','outside_capacity')#in litres
        self.parameters.vBS_ic=configuration.getArray('breeding_site','inside_capacity')#in litres
        self.parameters.vBS_od=configuration.getArray('breeding_site','outside_distribution')#distribution of BS outside # the sum of ouside and in-
        self.parameters.vBS_id=configuration.getArray('breeding_site','inside_distribution')#distribution of BS inside   #  side must be equal to 1
        self.parameters.vBS_os=configuration.getArray('breeding_site','outside_surface')#in cm^2
        self.parameters.n,self.parameters.m=len(self.parameters.vBS_od),len(self.parameters.vBS_id)
        self.parameters.ws_s=configuration.getFloat('weather','wind_shield')#wind shield in [0,1]
        #Cordoba
        self.parameters.location={'name':configuration.getString('location','name'),'station':configuration.getString('weather','station'),'zones':list(configuration.getString('location','zones'))}
        self.start_date=configuration.getDate('simulation','start_date')
        self.end_date=configuration.getDate('simulation','end_date')
        self.parameters.initial_condition=configuration.getArray('simulation','initial_condition')

        WEATHER_STATION_DATA_FILENAME='data/public/wunderground_'+self.parameters.location['station']+'.csv'
        self.parameters.weather=Weather(WEATHER_STATION_DATA_FILENAME,self.start_date,self.end_date)


    def getTimeRange(self):
        elapsed_days=(self.end_date - self.start_date).days
        return np.linspace(0, elapsed_days-1, elapsed_days * 20)

    def solveEquations(self,equations=equations.diff_eqs,method='odeint'):
        time_range = self.getTimeRange()
        initial_condition=self.parameters.initial_condition
        Y=None

        if(method=='odeint'):
            Y = spi.odeint(equations,initial_condition,time_range,hmax=1.0,args=(self.parameters,))#the ',' in (parameters,) is very important! '(parameters)' or tuple(parameters) doesn't work#TODO: this is because it calls aps out of it's domain.Find a better way.
        elif(method=='rk'):
            Y=rk.solve(equations,initial_condition,time_range,args=(self.parameters,))
        elif(method=='dopri'):
            Y=rk.scipy_solve(equations,initial_condition,time_range,'dopri',{'max_step':time_range[1]-time_range[0],'rtol':1e-3, 'atol':1e-6} )

        return time_range,initial_condition,Y
