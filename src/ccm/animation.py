#coding: utf-8
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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


def getCurves2D(fig,XYZ,X_M):
    ax = plt.axes()
    curve_M, = ax.plot(XYZ[:,0], XYZ[:,1])
    curve_X_M, = ax.plot(X_M[:,0], X_M[:,1])
    return curve_M,curve_X_M

def getCurves3D(fig,XYZ,X_M):
    ax=Axes3D(fig)
    curve_M, = ax.plot(XYZ[:,0], XYZ[:,1],XYZ[:,2])
    curve_X_M, = ax.plot(X_M[:,0], X_M[:,1],X_M[:,2])
    return curve_M,curve_X_M

mode='3D'

XYZ,time_range=solve()
#Animation code based on https://matplotlib.org/examples/animation/simple_anim.html
tau=5
X_M=np.array([ [XYZ[i-2*tau,0],XYZ[i,0],XYZ[i-tau,0]] for i in range(0,len(time_range)) ])
fig = plt.figure()
if(mode=='2D'):
    getCurves=getCurves2D
else:
    getCurves=getCurves3D

curve_M,curve_X_M=getCurves(fig,XYZ,X_M)

def animate2D(i):
    curve_M.set_data(XYZ[:i,0],XYZ[:i,1])  # update the data
    curve_X_M.set_data(X_M[:i,0],X_M[:i,1])  # update the data

def animate3D(i):
    animate2D(i)
    curve_M.set_3d_properties(XYZ[:i,2])
    curve_X_M.set_3d_properties(X_M[:i,2])

if(mode=='2D'):
    animate=animate2D
else:
    animate=animate3D

ani = animation.FuncAnimation(fig, animate, len(time_range), interval=1)

plt.show()
