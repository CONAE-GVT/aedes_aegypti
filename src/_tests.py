from otero_precipitation_wrapper import ModelWrapper as _Model
from config import Configuration#just importing these
from otero_precipitation import Model#two takes longer than solving _Model
import numpy as np


model=_Model("resources/otero_precipitation.cfg")
Y=np.array(model.solveEquations())
print(np.linalg.norm(Y))


configuration=Configuration('resources/otero_precipitation.cfg')
model=Model(configuration)
time_range,initial_condition,Y=model.solveEquations(method='rk' )
print(np.linalg.norm(Y))
