import pylab as pl
import numpy as np
from lorenz import solve
from scipy.stats import pearsonr

def exportToCsv(M,time_range):
    print('X, Y, Z')
    for i in range(0,len(time_range)):
        print('%s, %s, %s'%(M[i,0], M[i,1], M[i,2]))

#the d, and t this method returns seems to be ok, but it would be nice a second look
def nearestNeighbors(M_,time_range,E,tau):
    d=np.zeros( (len(time_range)- (E-1)*tau,E+1) )
    idx=np.zeros( (len(time_range)- (E-1)*tau,E+1),dtype=int )
    for i in range(0,len(time_range)-(E-1)*tau):
        distances_i=np.linalg.norm(M_-M_[i],axis=1)
        d[i,:]=np.sort(distances_i)[1:E+1 + 1]                        #[:E+1] -> [1: E+1 +1] to exclude distance to itself
        idx[i,:]=distances_i.argsort()[1:E+1 +1]#[:E+1] -> [1: E+1 +1] to exclude distance to itself

    return d,idx


def createWeights(d):
    w=np.zeros( (len(time_range)- (E-1)*tau,E+1) )
    for i in range(0,len(time_range)-(E-1)*tau):
        u=np.exp(-d[i,:]/d[i,0])
        N=np.sum(u)
        w[i,:]=u/N

    return w

def estimate(Y,M_,time_range,E,tau):
    d,idx=nearestNeighbors(M_,time_range,E,tau)
    w=createWeights(d)
    e=np.zeros( (len(time_range)- (E-1)*tau) )
    for i in range(0,len(time_range)-(E-1)*tau):
        e[i]=np.dot(w[i,:], Y[idx[i,:]+(E-1)*tau])

    return e#print(e-Y[(E-1)*tau:])

def C(Y,M_,time_range,E,tau):
    e=estimate(Y,M_,time_range,E,tau)
    rho,tail=pearsonr(Y[(E-1)*tau:],e)
    return rho**2

#x xmap y (y causes x)
#implementation follows: McCracken, James (2014). "Convergent cross-mapping and pairwise asymmetric inference".
if(__name__ == '__main__'):
    L_range=range(20,800,20)
    Z_crossmap_X=[]
    X_crossmap_Z=[]
    for L in L_range:
        M,time_range=solve(L)
        tau=1
        M_X=np.array([ [M[i,0],M[i-tau,0],M[i-2*tau,0]] for i in range(2*tau,len(time_range)) ])
        M_Y=np.array([ [M[i,1],M[i-tau,1],M[i-2*tau,1]] for i in range(2*tau,len(time_range)) ])
        M_Z=np.array([ [M[i,2],M[i-tau,2],M[i-2*tau,2]] for i in range(2*tau,len(time_range)) ])
        E=3
        #exportToCsv(M,time_range)
        #nearestNeighbors(M_X,time_range,E,tau)
        #createWeights(M_X,time_range,E,tau)
        #estimate(M[:,1],M_X,time_range,E,tau)
        X,Y,Z=M[:,0],M[:,1],M[:,2]
        Z_crossmap_X.append(C(X,M_Z,time_range,E,tau))
        X_crossmap_Z.append(C(Z,M_X,time_range,E,tau))
        print('.', end='', flush=True)

    pl.plot(L_range,Z_crossmap_X,label='Z xmap X')
    pl.plot(L_range,X_crossmap_Z,label='X xmap Z')
    pl.show()
