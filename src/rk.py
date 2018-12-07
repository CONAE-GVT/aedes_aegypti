import numpy as np
import pylab as pl
import math

def solve(_dYdt,Y0,time_range,args=(),steps=1):
    #main
    Y=np.zeros([len(time_range),len(Y0)])
    Y[0]=Y0#<---initial conditions
    def dYdt(Y,t,_args):
    	Y[Y<0]=0#this is to make rk work
    	return np.array(_dYdt(Y,t,*_args))#decorate the function to return an np array

    for i,t in enumerate(time_range[:-1]):
        h=time_range[i+1]-time_range[i]
        Y_j=Y[i]
        h_j=h/float(steps)
        for j in range(0,steps):
            #Runge-Kutta's terms
            K_n1=dYdt(Y_j,t,args)
            K_n2=dYdt(Y_j + (h_j/2.)*K_n1, t + h_j/2.,args)
            K_n3=dYdt(Y_j + (h_j/2.)*K_n2, t + h_j/2.,args)
            K_n4=dYdt(Y_j + h_j*K_n3, t + h_j,args)

            Y_j=Y_j+ (h_j/6.0)*(K_n1 + 2.0*K_n2 + 2.0*K_n3 + K_n4)
            t=t+h_j
        Y[i+1]=Y_j

    return Y

from scipy.integrate import ode
def scipy_solve(_dYdt,Y0,time_range,name,kwargs,args=()):
	def dYdt(t,Y): return np.array(_dYdt(Y,t,*args))#decorate the function to return an np array and swap args.
	r = ode(dYdt).set_integrator(name,**kwargs)
	r.set_initial_value(Y0, time_range[0])
	Y=np.zeros([len(time_range),len(Y0)])

	for i,t in enumerate(time_range[:-1]):
	    h=time_range[i+1]-time_range[i]
	    Y[i+1]=r.integrate(r.t+h)
	    if(not r.successful):
	        break

	return Y

def diff_eqs(Y,t):
	b = 0.25
	c = 5.0
	theta, omega = Y
	dydt = [omega, -b*omega - c*np.sin(theta)]
	return dydt

if (__name__ == '__main__'):
	Y0 = [np.pi - 0.1, 0.0]
	time_range=np.linspace(0, 10, 10*8)
	RES=solve(diff_eqs,Y0,time_range)

	pl.plot(RES[:,0])
	pl.plot(RES[:,1])
	pl.ylabel('y')
	pl.xlabel('t')
	pl.show()
