import kmeans
import numpy as np
import gdal
import sys
import os

FOLDER_IN='data/public/sensor/'
FOLDER_OUT='out/'


#gdal_translate -projwin 375873 6540024 398873 6512217 SPOT6.tiff SPOT6_cropped_.tiff
parameters={#TODO:check this coordinates and make them precise
    'west':'375873',
    'north':'6540024',
    'east':'398873',
    'south':'6512217'
 }

def gdalTranslate(parameters):
    #TODO: use gdal bindings
    cmd,out_filename=getGdalTranslateCmd(parameters)
    os.system(cmd)
    return out_filename


def getGdalTranslateCmd(parameters):
    in_filename=parameters['filename']
    projwin=parameters['west'] +" "+parameters['north'] +" " +parameters['east'] +" " +parameters['south']
    out_filename=FOLDER_OUT+'R.tif'
    cmd=(
    'gdal_translate -projwin {projwin} {in_filename} {out_filename}'
    ).format(projwin=projwin,in_filename=in_filename,out_filename=out_filename)
    return cmd,out_filename


def gdalMerge(parameters):
    #TODO: use gdal bindings
    cmd,out_filename=getGdalMergeCmd(parameters)
    os.system(cmd)
    return out_filename

def getGdalMergeCmd(parameters):
    out_filename=FOLDER_OUT+os.path.basename(parameters['filename'])
    in_filenames=parameters['in_filenames']
    cmd=(
    'rm {out_filename};gdal_merge.py -o {out_filename} -separate {in_filenames}'
    ).format(out_filename=out_filename,in_filenames=' '.join(in_filenames) )
    return cmd,out_filename

if(__name__=='__main__'):
    if(len(sys.argv)>2):#>python riskMap.py in/sentinel/B04.jp2 in/sentinel/B03.jp2 in/sentinel/B02.jp2
        parameters['filename']=FOLDER_OUT+'merged.tif'#out filename that will be created
        parameters['in_filenames']=sys.argv[1:]
        gdalMerge(parameters)
    else:
        parameters['filename']=sys.argv[1]

    print(kmeans.kmeans_classification(gdalTranslate(parameters),n_clusters=4))
