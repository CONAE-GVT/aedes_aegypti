#coding: utf-8
import sys
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
def lorenz(V,t):#parameters taken from http://node99.org/tutorials/ar/
    delta=10.
    beta=8/3.
    rho=28.
    x,y,z=V
    return [delta*(y-x),x*(rho-z)-y,x*y-beta*z]

def halvorsen(V,t):
    a=1.89
    x,y,z=V
    return [-a*x -4*y -4*z - y**2, -a*y - 4*z - 4*x -z**2,-a*z - 4*x - 4*y - x**2]

def diff_eqs(V,t):
    '''The main set of equations'''
    return lorenz(V,t)

def solve(L=100):
    initial_condition=[-8,8,27]
    time_range=np.linspace(0,L,L*100)
    return spi.odeint(diff_eqs,initial_condition,time_range),time_range


def getCurves(fig,M,M_X,M_Y):
    ax=Axes3D(fig)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    curve_M, = ax.plot(M[:,0], M[:,1],M[:,2],'-o',markevery=[-1],label='= (X(t), Y(t), Z(t))')
    curve_M_X, = ax.plot(M_X[:,0], M_X[:,1],M_X[:,2],'-o',markevery=[-1],label=r'= (X(t), X(t-$\tau$),X(t-2$\tau$))')
    curve_M_Y, = ax.plot(M_Y[:,0], M_Y[:,1],M_Y[:,2],'-o',markevery=[-1],label=r'= (Y(t), Y(t-$\tau$),Y(t-2$\tau$))')
    line,=ax.plot([],[])
    return curve_M,curve_M_X,curve_M_Y,line

def animate(i,M,M_X,M_Y,curve_M,curve_M_X,curve_M_Y,line):
    curve_M.set_data(M[:i,0],M[:i,1])  # update the data
    curve_M.set_3d_properties(M[:i,2])
    curve_M_X.set_data(M_X[:i,0],M_X[:i,1])  # update the data
    curve_M_X.set_3d_properties(M_X[:i,2])
    curve_M_Y.set_data(M_Y[:i,0],M_Y[:i,1])  # update the data
    curve_M_Y.set_3d_properties(M_Y[:i,2])
    line.set_data([ M_X[i-1,0],M_Y[i-1,0] ],[ M_X[i-1,1],M_Y[i-1,1] ])
    line.set_3d_properties([ M_X[i-1,2],M_Y[i-1,2] ])


if(__name__ == '__main__'):
    fig = plt.figure()
    M,time_range=solve()
    #Animation code based on https://matplotlib.org/examples/animation/simple_anim.html
    tau=2
    M_X=np.array([ [M[i,0],M[i-tau,0],M[i-2*tau,0]] for i in range(2*tau,len(time_range)) ])#+ np.array([20,0,0])
    M_Y=np.array([ [M[i,1],M[i-tau,1],M[i-2*tau,1]] for i in range(2*tau,len(time_range)) ])#+ np.array([-20,0,0])
    curve_M,curve_M_X,curve_M_Y,line=getCurves(fig,M,M_X,M_Y)
    ani = animation.FuncAnimation(fig, animate, range(0,len(time_range)), interval=100,fargs=(M,M_X,M_Y,curve_M,curve_M_X,curve_M_Y,line))
    if(len(sys.argv)>1 and  sys.argv[1]=='1'):
        for curve_M_ in [curve_M_X,curve_M_Y]:
            curve_M_.set_visible(False)
            curve_M_.set_label(None)
        line.set_visible(False)
    plt.legend(loc=0)
    plt.show()
