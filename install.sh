#!/bin/bash
sudo apt install -y python3-dev git python3-gdal python3-skimage python3-pip python3-grib python3-netcdf4
pip3 install scipy==0.19.1 numpy==1.13.3 matplotlib==2.1.1 similaritymeasures==0.3.0 scikit_image==0.13.1 plotly==2.2.3 folium==0.10.1 Pillow==7.0.0 scikit_learn==0.22.1 moviepy==0.2.3.5 --user
#cupy==5.1.0 pandas==0.22.0 optionals

#clone the project
git clone ssh://exequiel@pluton.hopto.org:8080/home/exequiel/repositories/aedes_aegypti.git

#**Installing pybind11**
git clone https://github.com/pybind/pybind11.git
sudo cp pybind11/include/pybind11 /usr/include/  -r

#**Installing Eigen**
git clone https://gitlab.com/libeigen/eigen.git
sudo cp eigen/Eigen /usr/include/  -r

#**c++ binding**
cd aedes_aegypti
g++ -std=c++14 -Wall -O3 -march=native -shared -fPIC -I/usr/include/python3.6m src/cpp/otero_precipitation_wrapper.cpp -o src/otero_precipitation_wrapper.so
