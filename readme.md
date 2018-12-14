**Dependencies**
>sudo apt install python3-scipy python3-grib python3-netcdf4 python3-matplotlib python3-gdal python3-xlrd python3-line-profiler
>pip install moviepy

**Run**
>python src/otero_precipitation.py
>python src/tests.py



**Git Setup**
On the remote server:
>mkdir aedes_aegypti.git
>cd aedes_aegypti.git/
>git --bare init --shared

On the local:
>git clone exequiel@zotac.hopto.org:/home/exequiel/repositories/aedes_aegypti.git
>touch readme.md
>git add readme.md
>git commit -m "initial commit"
>git push origin master
Total 18 (delta 0), reused 0 (delta 0)
To exequiel@zotac.hopto.org:/home/exequiel/repositories/aedes_aegypti.git
 * [new branch]      master -> master

*Not sure if I should on the local, create a repo and push it to the origin (instead of cloning, changing and pushing). like in https://feeding.cloud.geek.nz/posts/setting-up-centralied-git-repository/*


**To clone the project**
>git clone exequiel@zotac.hopto.org:/home/exequiel/repositories/aedes_aegypti.git
or
>git clone ssh://exequiel@zotac.hopto.org:8080/home/exequiel/repositories/aedes_aegypti.git

**To create a new branch**
>git checkout -b development
Commit if you have something
>git push -u origin development

**To delete a branch**
>git push origin --delete development
>git branch -d development

**To create a tag**
>git checkout master
>git tag -a v0.1 -m "Version: 0.1"
>git tag -l
>git push origin v0.1

**To merge**
First merge master -> development
>git checkout development
>git merge master development
>git push origin development

Now merge development -> master
>git checkout master
>git merge master development
>git push origin master

 *the branch order in the merge command is not important*

**Profiling**
>python -m cProfile -s cumtime src/tests.py  > p.txt
or
>kernprof -l src/tests.py spatial ; python -m line_profiler tests.py.lprof
source:https://github.com/rkern/line_profiler

**Usage**
To save:
>python src/tests.py save

To compare against previous results:
>python src/tests.py compare data/test/previous_results/2018-04-17__09_40_45.csv data/test/previous_results/2018-04-17__09_39_47.csv
or
>ls data/test/previous_results/*.csv |  python src/tests.py compare

**Cython**
>cd src
>python setup.py build_ext --inplace

**Installing Boost**
>wget https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.gz
>tar --gz -xf boost_1_67_0.tar.gz
>cd boost_1_67_0/
>sudo mkdir /home/exequiel/Downloads/libboost-python1.67/
>./bootstrap.sh --prefix=/home/exequiel/Downloads/libboost-python1.67/ -with-libraries=python
>./b2 install
>cd /home/exequiel/Downloads
inside libboost-python1.67, create a folder "usr" and move "lib" and "include" there.
Also create a folder "DEBIAN" with a file called "control" with Package:,Version:, Architecture:, etc. (you can download a deb package and use the control in there as a guide)
>dpkg-deb --build libboost-python1.67
>sudo dpkg -i libboost-python1.67
Source: http://www.king-foo.com/2011/11/creating-debianubuntu-deb-packages/


**c++ binding**
> g++ -Wall -std=c++11 -Wno-deprecated-declarations  -I /usr/local/boost_1_67_0/include/ -I/usr/include/python2.7 -fpic  src/_equationsmodule.cpp -shared  -L /usr/local/boost_1_67_0/lib/ -lboost_python27 -lboost_numpy27 -o src/_equations.so -O3
>export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/boost_1_67_0/lib
>python src/_equationsmodule_test.py

**To create an animation**
>ffmpeg -framerate 250 -i out/A_%04d.png  out/A.webm


**Engauge**
>pdfimages -j densityDependantlarvaeDeath.pdf /tmp/out

**Risk map**
>python src/risk_map.py data/public/sensor/sentinel/B02.jp2 data/public/sensor/sentinel/B03.jp2 data/public/sensor/sentinel/B04.jp2 data/public/sensor/sentinel/B08.jp2
band info: https://www.satimagingcorp.com/satellite-sensors/other-satellite-sensors/sentinel-2a/


**Notes**
constant ml -> constant survival (independent of initial larvaes)
