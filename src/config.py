import ConfigParser
import numpy as np
import datetime
class Configuration:
  def __init__(self, filename=None,dict_config=None):
    self.config_parser = ConfigParser.ConfigParser()
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

    #self.validate()

  def getFloat(self,section,option):
      return np.float(self.config_parser.get(section,option))

  def getArray(self,section,option):
      return np.fromstring(self.config_parser.get(section,option), dtype=float, sep=',')

  def getDate(self,section,option):
      return datetime.datetime.strptime(self.config_parser.get(section,option), '%Y-%m-%d').date()

  def getString(self,section,option):
      return self.config_parser.get(section,option)

  def validate(self):
    vBS_od=self.getArray('breeding_site','outside_distribution')
    vBS_id=self.getArray('breeding_site','inside_distribution')
    print(abs(1.-np.sum(vBS_od)+np.sum(vBS_id)))
    assert abs(1.-np.sum(vBS_od)+np.sum(vBS_id))<1e-10,'sum(vBS_id)+sum(vBS_id)!=1'

  def save(self,filename):
    self.config_parser.write(open(filename,'w'))
