import utils
import datetime
import numpy as np
import pylab as pl
from math import sin,cos,pi,exp,log
import fitter

def aFor(s,x0,P,n):
    return 2./P * sum([s(x) * (  sin(2.*pi*n*(x+1)/P)-sin(2.*pi*n*x/P) ) * 1./(2.*pi*n/P) for x in range(x0,P)] )

def bFor(s,x0,P,n):
    return 2./P * sum([s(x) * -( cos(2.*pi*n*(x+1)/P)-cos(2.*pi*n*x/P) ) * 1./(2.*pi*n/P) for x in range(x0,P)] )

def fourier(real_values,N):
    x0=0#TODO:not sure about this
    P=len(real_values)#period
    s=lambda t: real_values[t]#step function for real_values
    a0=2./P* sum([s(x) for x in range(x0,P)])
    a=np.array([a0]+[aFor(s,x0,P,n) for n in range(1,N)])
    b=np.array([0]+[bFor(s,x0,P,n) for n in range(1,N)])
    n=np.array(range(1,N))
    return lambda x:a[0]/2. + np.sum(a[1:N]* np.cos(2.*pi*n*x/P) + b[1:N]* np.sin(2.*pi*n*x/P) )

if (__name__ == '__main__'):
    real_values=[10,15,13,16,15,15,14,18,19,10]
    s_5=fourier(real_values,5)
    #pl.plot([s(x) for x in range(x0,P)], '-y')    
    #pl.plot([s_5(x) for x in range(x0,P)], '-b')



    start_date=datetime.date(2014, 7, 1)
    end_date=datetime.date(2017, 6, 1)
    WEATHER_STATION_DATA_FILENAME='data/wunderground_SACO.csv'
    temperatures= utils.getAverageTemperaturesFromCsv(WEATHER_STATION_DATA_FILENAME,start_date,end_date)
    T=fourier(temperatures,50)
    pl.plot(temperatures,'-y', label='Data')#.y', label='Data',markersize=3    
    a,b,c=fitter.getOptimalParameters(temperatures)
    pl.plot([fitter.T(a,b,c,t) for t in range(0,len(temperatures))], label='Usual T(t)')
    pl.plot([T(t) for t in range(0,len(temperatures))],'-m', label='Fourier')
    pl.legend(loc=0)
    
    #wind_speeds=utils.getMeanWindSpeedFromCsv(WEATHER_STATION_DATA_FILENAME,start_date,end_date)    
    #ws=fourier(wind_speeds,50)
    #pl.plot([step_ws(t) for t in range(0,len(wind_speeds))], '-y')
    #pl.plot([ws(t) for t in range(0,len(wind_speeds))], '-b')

    #interval=np.linspace(1,0,200)
    #pl.plot(interval,[(1-exp(z))/(1-exp(1)) for z in interval])
    #pl.plot(interval,[log(1+(exp(1)-1)*z**0.25) for z in interval])

    pl.show()
    print("||old_fourierT-current_fourierT|| = "),
    old_fourierT=np.load('old_fourierT.npy')
    print(np.linalg.norm(old_fourierT-[T(t) for t in range(0,len(temperatures))]))
    #np.save('backup/RES_%s.npy'%(datetime.datetime.now().isoformat()),RES)

