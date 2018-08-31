#coding: utf-8
import sys
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


def getCurves2D(fig,M,M_X):
    ax = plt.axes()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    curve_M, = ax.plot(M[:,0], M[:,1])
    curve_M_X, = ax.plot(M_X[:,0], M_X[:,1])
    return curve_M,curve_M_X

def getCurves3D(fig,M,M_X):
    ax=Axes3D(fig)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    curve_M, = ax.plot(M[:,0], M[:,1],M[:,2],'-o',markevery=[-1],label='= (X(t), Y(t), Z(t))')
    curve_M_X, = ax.plot(M_X[:,0], M_X[:,1],M_X[:,2],'-o',markevery=[-1],label=r'= (X(t-2$\tau$), X(t), X(t-$\tau$))')
    return curve_M,curve_M_X

def animate2D(i,M,M_X,curve_M,curve_M_X):
    curve_M.set_data(M[:i,0],M[:i,1])  # update the data
    curve_M_X.set_data(M_X[:i,0],M_X[:i,1])  # update the data

def animate3D(i,M,M_X,curve_M,curve_M_X):
    animate2D(i,M,M_X,curve_M,curve_M_X)
    curve_M.set_3d_properties(M[:i,2])
    curve_M_X.set_3d_properties(M_X[:i,2])


if(__name__ == '__main__'):
    mode='3D'
    if(mode=='2D'):
        getCurves=getCurves2D
        animate=animate2D
    else:
        getCurves=getCurves3D
        animate=animate3D

    fig = plt.figure()
    M,time_range=solve()
    #Animation code based on https://matplotlib.org/examples/animation/simple_anim.html
    tau=5
    M_X=np.array([ [M[i-2*tau,0],M[i,0],M[i-tau,0]] for i in range(2*tau,len(time_range)) ])
    curve_M,curve_M_X=getCurves(fig,M,M_X)
    ani = animation.FuncAnimation(fig, animate, range(3000,len(time_range)), interval=1,fargs=(M,M_X,curve_M,curve_M_X))
    if(len(sys.argv)>1 and  sys.argv[1]=='1'):
        curve_M_X.set_visible(False)
        curve_M_X.set_label(None)
    plt.legend(loc=0)
    plt.show()
