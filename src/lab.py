import numpy as np
import utils
import rk
import plotly.graph_objs as go

def runTest():
    W0=3
    x=np.linspace(0,20,20*5)
    y=np.exp(-(x-W0)**2)/np.sqrt(np.pi)
    print((y* 20./len(y)).sum())
    utils.showPlot([go.Scatter(x=x,y=y, name='')],title='')

def dW(t):
    return np.sin(2*np.pi*t/2)
    if(6.5<t<7.0):
        return 20e-1
    else:
        return -3e-1

def dE(E,W,t):
    ovsp_t=0.5
    h=2
    return ovsp_t*np.exp(-(2-W)**2)/np.sqrt(np.pi)#t-W makes no sense! it's h-w integrated over h. now it's just how much eggs accumulated at level h

def diff_eqs(Y,t,h,parameters):
    '''The main set of equations'''
    # T_t=parameters.weather.T(t)
    # p_t=parameters.weather.p(t)
    # RH_t=parameters.weather.RH(t)
    # elr,lpr,par,ovr1,ovr2=vR_D(T_t)

    #dY=np.empty(Y.shape)
    E,W=Y
    return np.array([dE(E,W,t),dW(t)])


def solve():
    initial_condition=[5,2]
    time_range=np.linspace(0,10,40)
    Y=rk.solve(diff_eqs,initial_condition,time_range,args=([],),steps=20)
    utils.showPlot([go.Scatter(x=time_range,y=Y[:,0], name='E'),
                    go.Scatter(x=time_range,y=Y[:,1], name='W')],title='')

if(__name__ == '__main__'):
    solve()
