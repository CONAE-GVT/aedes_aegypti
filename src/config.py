from configparser import ConfigParser
import numpy as np
import datetime
class Configuration:
  def __init__(self, filename=None,dict_config=None):
    self.config_parser = ConfigParser()
    if(filename):
        self.config_parser.readfp(open(filename))
    if(dict_config):
        for section in dict_config:
            for option in dict_config[section]:
                values=dict_config[section][option]
                if(isinstance(values,str)):
                    string_value=values
                else:#not string, try to convert it
                    try:#it's iterable
                        string_value=','.join([str(value) for value in values ])
                    except TypeError:
                        string_value=str(values)
                self.config_parser.set(section,option,string_value)

    self.validate()

  def getFloat(self,section,option):
      return np.float(self.config_parser.get(section,option))

  def getArray(self,section,option):
      return np.fromstring(self.config_parser.get(section,option), dtype=float, sep=',')

  def getDate(self,section,option):
      return datetime.datetime.strptime(self.config_parser.get(section,option), '%Y-%m-%d').date()

  def getString(self,section,option):
      return self.config_parser.get(section,option)

  def validate(self):
    vBS_d=self.getArray('breeding_site','distribution')
    vBS_s=self.getArray('breeding_site','surface')
    vBS_h=self.getArray('breeding_site','height')
    vBS_mf=self.getArray('breeding_site','manually_filled')
    initial_condition=self.getArray('simulation','initial_condition')
    alpha0=self.getArray('biology','alpha0')
    assert np.all(vBS_d>0),'vBS_d cannot have a zero in it'#not allowed anymore
    assert len(vBS_d) == len(vBS_h) == len(vBS_s) == len(vBS_mf),'vBS_d, vBS_h, vBS_s and vBS_mf must have the same dimension!'
    n=len(vBS_d)
    assert len(alpha0)== n, 'dim(alpha0)!=%s'%(n)
    assert abs(1.-np.sum(vBS_d))<1e-10,'sum(vBS_d)=%s!=1'%(np.sum(vBS_d))
    assert len(initial_condition)==3*n+1+1+n#(vE+vL+vP) +A1 +A2 +vW


  def save(self,filename):
    self.config_parser.write(open(filename,'w'))
