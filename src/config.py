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
            #self.config_parser.add_section(section)
            for option in dict_config[section]:
                values=dict_config[section][option]
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
    vBS_oc=self.getArray('breeding_site','outside_capacity')
    vBS_ic=self.getArray('breeding_site','inside_capacity')
    vBS_od=self.getArray('breeding_site','outside_distribution')
    vBS_id=self.getArray('breeding_site','inside_distribution')
    vBS_os=self.getArray('breeding_site','outside_surface')
    initial_condition=self.getArray('simulation','initial_condition')
    alpha0=self.getArray('biology','alpha0')
    assert np.all(vBS_id>0),'vBS_id cannot have a zero in it'#not allowed anymore, to achieve the same,
    assert np.all(vBS_od>0),'vBS_od cannot have a zero in it'#just don't put the container (empty arrays are allowed)
    assert len(vBS_od) == len(vBS_oc) == len(vBS_os),'vBS_od,vBS_oc and vBS_os must have the same dimension!'
    assert len(vBS_id) == len(vBS_ic),'vBS_id and vBS_ic must have the same dimension!'
    n,m=len(vBS_od),len(vBS_id)
    assert len(alpha0)== n+m, 'dim(alpha0)!=%s'%(n+m)
    assert len(initial_condition)==3*(n+m)+1+1+n#(vE+vL+vP) +A1 +A2 +vW
    assert abs(1.-np.sum(vBS_od)-np.sum(vBS_id))<1e-10,'sum(vBS_id)+sum(vBS_id)=%s!=1'%(np.sum(vBS_od)+np.sum(vBS_id))

  def save(self,filename):
    self.config_parser.write(open(filename,'w'))
