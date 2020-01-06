import os
import tempfile
import datetime
import numpy as np
from otero_precipitation import Model as pyModel
from config import Configuration

try:
    from otero_precipitation_wrapper import ModelWrapper as _Model
except ImportError:
    os.system('g++ -std=c++17 -O3 -march=native -Wall -I/usr/include/python3.6m/ -fpic src/cpp/otero_precipitation_wrapper.cpp -shared -lboost_python-py36 -o src/otero_precipitation_wrapper.so')
    from otero_precipitation_wrapper import ModelWrapper as _Model

class Model:
    config_filename=tempfile.NamedTemporaryFile(suffix='.cfg').name

    def __init__(self, configuration=Configuration('resources/otero_precipitation.cfg')):
        with open(self.config_filename, 'w') as configfile:
            configuration.config_parser.write(configfile)

        self._model=_Model(self.config_filename)
        self.start_date=datetime.datetime.strptime(self._model.start_date,'%Y-%m-%d').date()
        self.end_date=datetime.datetime.strptime(self._model.end_date,'%Y-%m-%d').date()
        self.time_range=np.array(self._model.time_range)
        self.warnings=['warnings not implemented']
        #TODO:we are cheating here!
        python_model=pyModel(configuration)
        self.parameters=python_model.parameters

    def solveEquations(self):
        self._model.solveEquations()
        self.Y=np.array(self._model.Y)
        return self.time_range,self.Y
