from scipy import optimize
import numpy as np
import pylab as pl
import datetime
import utils
import math

def f(a,b,c,z):
    #print(a- b*z/256.)
    #return math.log(a- b*z/256.)+c
    #return a*gamma.pdf(z,b) +c
    return 1-a*math.exp(b*z)/(c+math.exp(b*z))#inspired by x'=x*(1-x/k), and equvalent to b/(b+exp(a*z)) using a=1 which makes the function positive for all z

def f_error(x,domain,real_values):
    error= real_values - [ f(x[0],x[1],x[2],z) for z in domain]
    return np.dot(error,error)#for leastsq return just error

def getOptimalParameters(domain,real_values):
    x0 = np.array([0.0,0.0,0.])#initial a,b and c

    constraints = ({'type': 'ineq', 'fun': lambda x:  np.min(x[0]- x[1]*domain/256.)-1e-10})
    bounds=((0,2),(-1,1),(0,1))
    #res,cov=optimize.leastsq(f_error,x0,(domain,real_values))
    res=optimize.minimize(f_error,x0,(domain,real_values),method='SLSQP')#,constraints=constraints
    #res=optimize.differential_evolution(f_error,bounds=bounds,args=(domain,real_values))
    return res

if(__name__ == '__main__'):
    #define values to fit
    domain=     np.array([4,8   ,16  ,32  ,64    ,128  ,256])
    mortality=1-np.array([1,1   ,1   ,1   ,0.970 ,0.450,0.065])#mortality=1-survival in [0,1]
    pupation=np.array([5.45,5.78,5.41,7.14,11.75 ,10.39,7.42])#pupation time in days
    #find optimal
    res=getOptimalParameters(domain,mortality)
    a,b,c=res.x
    #print([a,b,c])
    print(res)


     #Ploting
    pl.plot(range(0,domain.max()),[f(a,b,c,x) for x in range(0,domain.max())],label='f')
    pl.plot(domain,mortality, '^y',label='')
    pl.xlabel('')
    pl.ylabel('')
    pl.legend(loc=0)

    pl.show()
