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

def getValues(location):
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
    values={}
    dates=[datetime.date(2018,12,1),datetime.date(2019,1,1),datetime.date(2019,2,1),datetime.date(2019,3,1),datetime.date(2019,4,1)]
    start_date=configuration.getDate('simulation','start_date')
    for i in range(len(dates)-1):
        date,next_date=dates[i],dates[i+1]
        start,end=(date-start_date).days,(next_date-start_date).days
        values[date.month]=np.max(lwO[indexOf(start):indexOf(end)])
    return values


def generateCSV():
    #return DATA_FOLDER+'arg_cities.csv'
    config_parser = ConfigParser()
    config_parser.read('resources/get_weather.cfg')
    params=[]
    start_date,end_date=datetime.date(2015,7,15),datetime.date(2019,5,1)
    output='location,month,value,lat,lon\n'
    for location in config_parser.sections():
        lat=float(config_parser.get(location,'lat'))
        lon=float(config_parser.get(location,'lon'))
        values=getValues(location)
        for month in values:
            output+='%s, %s, %s, %s, %s\n'%(location,month,values[month],lat,lon)
    filename=DATA_FOLDER+'arg_cities.csv'
    open(filename,'w').write(output)
    return filename

def plotMap():
    df = pd.read_csv(generateCSV())
    df.head()

    cases = []
    colors = ['rgb(239,243,255)','rgb(189,215,231)','rgb(107,174,214)','rgb(33,113,181)']
    months = {12:'Dec',1:'Jan',2:'Feb',3:'Mar'}

    for i in [12,1,2,3][::-1]:
        cases.append(go.Scattergeo(
                lon = df[ df['month'] == i ]['lon'], #-(max(range(6,10))-i),
                lat = df[ df['month'] == i ]['lat'],
                text = df[ df['month'] == i ]['value'],
                name = months[i],
                marker = go.scattergeo.Marker(
                    size = df[ df['month'] == i ]['value']/25,
                    color = colors[i%4],
                    line = go.scattergeo.marker.Line(width = 0)
                )
            )
        )

    cases[0]['text'] = df[ df['month'] == 1 ]['value'].map('{:.0f}'.format).astype(str)+' '+\
        df[ df['month'] == 1 ]['location']
    cases[0]['mode'] = 'markers+text'
    cases[0]['textposition'] = 'bottom center'

    inset = [
        go.Choropleth(
            locationmode = 'country names',
            locations = df[ df['month'] == 1 ]['location'],
            z = df[ df['month'] == 1 ]['value'],
            text = df[ df['month'] == 9 ]['location'],
            colorscale = [[0,'rgb(0, 0, 0)'],[1,'rgb(0, 0, 0)']],
            autocolorscale = False,
            showscale = False,
            geo = 'geo2'
        ),
        go.Scattergeo(
            lon = [21.0936],
            lat = [7.1881],
            text = ['Argentina'],
            mode = 'text',
            showlegend = False,
            geo = 'geo2'
        )
    ]

    layout = go.Layout(
        title = go.layout.Title(
            text = 'Ebola cases reported by month in West Africa 2014<br> \
    Source: <a href="https://data.hdx.rwlabs.org/dataset/rowca-ebola-cases">\
    HDX</a>'),
        geo = go.layout.Geo(
            resolution = 50,
            scope = 'south america',
            showframe = False,
            showcoastlines = True,
            showland = True,
            landcolor = "rgb(229, 229, 229)",
            countrycolor = "rgb(255, 255, 255)" ,
            coastlinecolor = "rgb(255, 255, 255)",
            projection = go.layout.geo.Projection(
                type = 'mercator'
            )
        ),
        legend = go.layout.Legend(
               traceorder = 'reversed'
        )
    )

    fig = go.Figure(layout=layout, data=cases)
    ply.plot(fig, filename='2017-2018 maximum weekly oviposition')

if(__name__ == '__main__'):
    plotMap()
