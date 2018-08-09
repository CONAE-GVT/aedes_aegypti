#coding: utf-8
import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation



def diff_eqs(V,t):
    '''The main set of equations'''
    delta=10.
    beta=8/3.
    rho=28.
    x,y,z=V
    return [delta*(y-x),x*(rho-z)-y,x*y-beta*z]

def solve():
    initial_condition=[1.,0,0]
    time_range=np.linspace(0,100,100**2)
    return spi.odeint(diff_eqs,initial_condition,time_range),time_range


XYZ,time_range=solve()
#Animation code based on https://matplotlib.org/examples/animation/simple_anim.html
fig, ax = plt.subplots()
curve_M, = ax.plot(XYZ[:,0], XYZ[:,1])
curve_X_M, = ax.plot(XYZ[0,0], XYZ[0,0],'o')

def animate(i):
    curve_M.set_data(XYZ[:i,0],XYZ[:i,1])  # update the data
    tau=5
    curve_X_M.set_data(XYZ[i-2*tau,0],XYZ[i,0])  # update the data
    return curve_M

ani = animation.FuncAnimation(fig, animate, len(time_range), interval=1)

plt.show()
