import gdal
import numpy as np
from sklearn.cluster import MiniBatchKMeans

def createOutputImage(out_filename, in_dataset):
    driver = gdal.GetDriverByName( "GTiff" )
    # Create an output file of the same size as the inputted image but with only 1 output image band.
    newDataset = driver.Create(out_filename, in_dataset.RasterXSize,in_dataset.RasterYSize, 1, gdal.GDT_Float32)
    # Define the spatial information for the new image.
    newDataset.SetGeoTransform(in_dataset.GetGeoTransform())
    newDataset.SetProjection(in_dataset.GetProjection())
    return newDataset

def kmeans_classification(filename,n_clusters=2):
    #open input data set
    in_dataset=gdal.Open(filename, gdal.GA_ReadOnly)
    raster=in_dataset.ReadAsArray()
    #save input data set as npy
    np.save('out/R.npy',raster)

    #create output dataset (classified)
    out_dataset=createOutputImage('out/C.tif',in_dataset)
    '''
    The original is a 3-dimensional matrix [band][x][y]
    for example, a matrix of [4,460,380]
    [4,460,380]---reshape--->[4,460*380]---transpose-->[[460*380,4]]
    '''
    flat_raster=raster.reshape(raster.shape[0],raster.shape[1]*raster.shape[2]).transpose()

    kmeans=MiniBatchKMeans(n_clusters=n_clusters).fit(flat_raster)
    code=kmeans.predict(flat_raster.astype(float))#as float to avoid a warning.
    code_image=code.reshape(raster.shape[1],raster.shape[2])#this are class numbers, integers from [0,n_clusters]
    #save the classified raster as tif and npy
    out_dataset.GetRasterBand(1).WriteArray(code_image)#TODO:maybe we should use a color map like envi. >gdalinfo out/SPOT6.kmeans.envi.tiff
    np.save('out/C.npy',code_image)

    #
    print('Cluster centers:\n'+str(kmeans.cluster_centers_))
    print('Code Image(image of labeled pixels):\n'+str(code_image))
    print('Inertia:'+str(kmeans.inertia_))
    print('Min: '+ str(min(code_image.flatten()))+ ' Max: '+str(max(code_image.flatten())))
    #

    return 'out/C.tif'
