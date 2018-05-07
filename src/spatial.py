from otero_precipitation import Model
from config import Configuration
from _equations import diff_eqs
import numpy as np
import pylab as pl
import matplotlib
import datetime
import sys


WIDTH=20
HEIGHT=20
FRAMES=10
def getConfig(configuration,i,j):
    vBS_od=configuration.getArray('breeding_site','outside_distribution')
    vBS_id=configuration.getArray('breeding_site','inside_distribution')

    vBS_od=np.array([1. for x in vBS_od])#vBS_od/np.sum(vBS_od)
    vBS_id=np.array([1. for x in vBS_id])#vBS_id/np.sum(vBS_id)

    BS_o=float(i+j+1)/float(WIDTH+HEIGHT)
    vBS_od=BS_o*vBS_od
    vBS_id=(1-BS_o)*vBS_id
    assert 1-(np.sum(vBS_od)+np.sum(vBS_id))<1e-6, 'vBS_od+vBS_id=%s'%(np.sum(vBS_od)+np.sum(vBS_id))

    configuration.config_parser.set('breeding_site','outside_distribution',','.join([str(value) for value in vBS_od ]))
    configuration.config_parser.set('breeding_site','inside_distribution',','.join([str(value) for value in vBS_id ]))
    return configuration

if(__name__ == '__main__'):
    config=Configuration('resources/otero_precipitation.cfg',{
        'simulation':{
            'start_date':datetime.date(2017,10,1),
            'end_date':datetime.date(2018,4,5)
        }
    })

    matrix=None
    for i in range(0,HEIGHT):
        for j in range(0,WIDTH):
            model=Model(getConfig(config,i,j))
            time_range,initial_condition,Y=model.solveEquations(equations=diff_eqs,method='rk')
            if matrix is None:
                matrix=np.zeros((WIDTH,HEIGHT,len(time_range)))

            ADULT1,ADULT2=model.parameters.ADULT1,model.parameters.ADULT2
            matrix[i,j,:]=Y[:,ADULT1]+Y[:,ADULT2]

            sys.stdout.write("Progress: %d%%   \r" % (float(i*HEIGHT+j)/float(WIDTH*HEIGHT) * 100.) )
            sys.stdout.flush()

    '''
    #show several plots
    normalizer=matplotlib.colors.Normalize(vmin=matrix.min(), vmax=matrix.max(), clip=False)
    for a in range(0,FRAMES):
        size=len(time_range)
        n=int((a/float(FRAMES))*float(size))
        pl.figure()
        pl.imshow(matrix[:,:,n],norm=normalizer,cmap='gray')
    pl.show()
    '''

    raw_input('Press enter')
    #https://matplotlib.org/gallery/animation/animation_demo.html#sphx-glr-gallery-animation-animation-demo-py
    matrix[:,:,i]=matrix[:,:,i]/matrix.max()
    fig, ax = pl.subplots()
    for i in range(0,len(time_range)):
        ax.cla()
        ax.imshow(matrix[:,:,i],cmap='gray',interpolation="nearest")
        ax.set_title("frame {}".format(i))
        pl.pause(0.0001)
