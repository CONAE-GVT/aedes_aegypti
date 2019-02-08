from otero_precipitation_wrapper import ModelWrapper as _Model
from config import Configuration#just importing these
from otero_precipitation import Model#two takes longer than solving _Model
import numpy as np

if(__name__ == '__main__'):
    model=_Model('resources/otero_precipitation.cfg')
    Y1=np.array(model.solveEquations())
    print(np.linalg.norm(Y1),Y1.shape)


    configuration=Configuration('resources/otero_precipitation.cfg')
    model=Model(configuration)
    time_range,initial_condition,Y2=model.solveEquations(method='rk' )
    print(np.linalg.norm(Y2),Y2.shape)

    print(Y1[:,0])
    print(Y2[:,0])
