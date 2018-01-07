from scipy import optimize
import numpy as np
import pylab as pl
import datetime
import utils
import math



#real_temps =utils.getAverageTemperaturesFromCsv('data/noaa_big.csv',datetime.date(2001, 07, 1),datetime.date(2017, 6, 1))#data for 2003-09-23 was missing, so min+max/2 was put in it-s place
#real_temps =textReader.getAverageTemperaturesFromTxts('data/temperatures',datetime.date(2016, 07, 1),datetime.date(2017, 6, 1))#np.array([26.0,30.0,23.0])

def T(a,b,c,t,noise=0):#K, beginning on the first of july. Noise is just used to test the fitter. For the actual use of the function, the noise is always 0
    #a=18.0#C
    #b=6.7#C
    #c=9.2
    return a + b * math.cos( (2.0* math.pi * t)/365.25 + c) + 273.15 +noise#To Kelvin

#Just to test
#vT=np.vectorize(T)
#real_temps = vT(20.0,10.0,15.0,np.arange(0, 365*10, 1.0),np.random.normal(0,2,365*10))#np.array([26.0,30.0,23.0])
def T_error(x,real_temps):
    return [ real_temps[i]-T(x[0],x[1],x[2],i) if real_temps[i] else 0 for i in range(0,len(real_temps))]

def getOptimalParameters(real_temps):
    x0 = np.array([18.0,6.7,9.2])#initial a,b and c
    res,cov=optimize.leastsq(T_error,x0,(real_temps))
    return res

if(__name__ == '__main__'):
    real_temps=utils.getAverageTemperaturesFromCsv('data/noaa_big.csv',datetime.date(2001, 7, 1),datetime.date(2017, 6, 1))
    a,b,c=getOptimalParameters(real_temps)
    print([a,b,c])

    #Ploting
    time_range = np.arange(0, len(real_temps), 1.0)
    vT=np.vectorize(T)
    pl.subplot(111)
    pl.plot(time_range,vT(a,b,c,time_range), '-m', label='T')
    pl.xlabel('Time(in days starting in July)')
    pl.ylabel('')
    pl.legend(loc=0)



    pl.subplot(111)
    pl.plot(time_range,real_temps, '^y', label='Real Temperatures')
    pl.legend(loc=0)


    pl.show()
