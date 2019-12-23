**Dependencies**
>sudo apt install python3-scipy python3-grib python3-netcdf4 python3-matplotlib python3-gdal python3-line-profiler nvidia-cuda-toolkit python3-skimage python3-plotly ffmpeg python3-msgpack python3-pandas
>sudo -H pip install moviepy cupy
>pip install plotly --upgrade --user
**Run**
>python src/otero_precipitation.py
>python src/tests.py



**Git Setup**
On the remote server:
>mkdir aedes_aegypti.git
>cd aedes_aegypti.git/
>git --bare init --shared

On the local:
>git clone exequiel@pluton.hopto.org:/home/exequiel/repositories/aedes_aegypti.git
>touch readme.md
>git add readme.md
>git commit -m "initial commit"
>git push origin master
Total 18 (delta 0), reused 0 (delta 0)
To exequiel@pluton.hopto.org:/home/exequiel/repositories/aedes_aegypti.git
 * [new branch]      master -> master

*Not sure if I should on the local, create a repo and push it to the origin (instead of cloning, changing and pushing). like in https://feeding.cloud.geek.nz/posts/setting-up-centralied-git-repository/*


**To clone the project**
>git clone exequiel@pluton.hopto.org:/home/exequiel/repositories/aedes_aegypti.git
or
>git clone ssh://exequiel@pluton.hopto.org:8080/home/exequiel/repositories/aedes_aegypti.git

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
>sudo apt install libboost-python1.65-dev

**c++ binding**
>g++ -O3 -Wall -I/usr/include/python3.6m/   -fpic  src/cpp/otero_precipitation_wrapper.cpp -shared  -lboost_python-py36 -lpython3.6m -o src/otero_precipitation_wrapper.so

**To create an animation**
>ffmpeg -framerate 250 -i out/A_%04d.png  out/A.webm


**Engauge**
>pdfimages -j densityDependantlarvaeDeath.pdf /tmp/out

**Risk map**
>python src/risk_map.py data/public/sensor/sentinel/B04.jp2 data/public/sensor/sentinel/B03.jp2 data/public/sensor/sentinel/B02.jp2 data/public/sensor/sentinel/B08.jp2
band info: https://www.satimagingcorp.com/satellite-sensors/other-satellite-sensors/sentinel-2a/



**Web**
>sudo apt install apache2 libapache2-mod-python
>sudo a2enmod cgid
>sudo service apache2 restart
>sudo vim /etc/apache2/conf-available/cgi-enabled.conf
# create new
# processes .cgi and .py as CGI scripts
<Directory "/var/www/html/app/cgi-bin">
   Options +ExecCGI
   AddHandler cgi-script .cgi .py
</Directory>

>sudo ln -s /home/exequiel/Desktop/models/programs/aedes_aegypti/src/app /var/www/html/
>sudo a2enconf cgi-enabled
>sudo service apache2 restart
>chmod 705 /home/exequiel/Desktop/models/programs/aedes_aegypti/src/app/cgi-bin/iface.py
#optional:
>sudo vim /etc/apache2/conf-available/charset.conf #and uncomment the utf-8
>sudo mkdir /var/www/.imageio ; sudo chmod 777 /var/www/.imageio -R #to allow moviepy download ffmpeglib
#Can now access to http://localhost/app/cgi-bin/iface.py

**Notes**
constant ml -> constant survival (independent of initial larvaes)


#doc
IEEE template
>sudo apt install texlive-publishers

QGIS:
Open arg.jp2
add layer> add delimited text. (select cities.csv)
right click cities layer and click properties:
  on the left select "Labels":
    In Text>"label with" select city
    in Buffer > draw buffer
source: http://transnationalhistory.net/mapping/tutorials/pointvectorlayers/


#AT plots
python src/tests.py 18;python src/tests.py 8;python src/tests.py 8;python src/tests.py 11;python src/tests.py 12;python src/tests.py 14;
