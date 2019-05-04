from configparser import ConfigParser
from otero_precipitation import Model
from equations import diff_eqs
from config import Configuration
import utils
import numpy as np
import plotly.graph_objs as go
import plotly.offline as ply
import pandas as pd
import datetime


DATA_FOLDER='data/public/'

def getValue(location):
    #run simulation for location
    mf=0.
    h=10.
    configuration=Configuration('resources/2c.cfg')
    configuration.config_parser.set('location','name',location+'.full')
    configuration.config_parser.set('simulation','end_date',str(datetime.date.today()))
    n=len(configuration.getArray('breeding_site','height'))
    configuration.config_parser.set('breeding_site','height',','.join([str(h)]*n))
    configuration.config_parser.set('breeding_site','manually_filled',','.join([str(mf)]+[str(0)]*(n-1)))
    model=Model(configuration)
    time_range,initial_condition,Y=model.solveEquations(equations=utils.OEquations(model,diff_eqs),method='rk')
    #calculate oviposition
    indexOf=lambda t: (np.abs(time_range-t)).argmin()
    OVIPOSITION=model.parameters.OVIPOSITION
    BS_a=model.parameters.BS_a
    O=Y[:,OVIPOSITION]
    lwO=np.sum([Y[indexOf(t),OVIPOSITION]-Y[indexOf(t-7),OVIPOSITION] for t in time_range],axis=1)/BS_a#calculate the difference,sum up all levels, and divide by amount of containers
    return np.max(lwO)


def generateCSV():
    #return DATA_FOLDER+'arg_cities.csv'
    config_parser = ConfigParser()
    config_parser.read('resources/get_weather.cfg')
    params=[]
    start_date,end_date=datetime.date(2015,7,15),datetime.date(2019,5,1)
    output='name,value,lat,lon\n'
    for location in config_parser.sections():
        lat=float(config_parser.get(location,'lat'))
        lon=float(config_parser.get(location,'lon'))
        value=getValue(location)
        output+='%s, %s, %s, %s\n'%(location,value,lat,lon)
    filename=DATA_FOLDER+'arg_cities.csv'
    open(filename,'w').write(output)
    return filename

def plotMap():
    df = pd.read_csv(generateCSV())
    df.head()

    df['text'] = df['name'] + '<br>Max weekly oviposition ' + df['value'].astype(str)
    limits = [(0,2),(3,10),(11,20),(21,50),(50,3000)]
    colors = ["rgb(0,116,217)","rgb(255,65,54)","rgb(133,20,75)","rgb(255,133,27)","lightgrey"]
    cities = []


    for i in range(len(limits)):
        lim = limits[i]
        df_sub = df[lim[0]:lim[1]]
        city = go.Scattergeo(
            locationmode = 'country names',
            lon = df_sub['lon'],
            lat = df_sub['lat'],
            text = df_sub['text'],
            marker = go.scattergeo.Marker(
                size = df_sub['value'],
                color = colors[i],
                line = go.scattergeo.marker.Line(
                    width=0.5, color='rgb(40,40,40)'
                ),
                sizemode = 'area'
            ),
            name = '{0} - {1}'.format(lim[0],lim[1]) )
        cities.append(city)

    layout = go.Layout(
            title = go.layout.Title(
                text = '2018-2019 Argentina'
            ),
            showlegend = True,
            geo = go.layout.Geo(
                scope = 'south america',
                projection = go.layout.geo.Projection(
                    type='mercator'
                ),
                showland = True,
                landcolor = 'rgb(217, 217, 217)',
                subunitwidth=1,
                countrywidth=1,
                subunitcolor="rgb(255, 255, 255)",
                countrycolor="rgb(255, 255, 255)"
            )
        )

    fig = go.Figure(data=cities, layout=layout)
    ply.plot(fig, filename='d3-bubble-map-populations')


if(__name__ == '__main__'):
    plotMap()
