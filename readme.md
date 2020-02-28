**Dependencies**
>sudo apt install python3-scipy python3-grib python3-netcdf4 python3-matplotlib python3-gdal python3-line-profiler nvidia-cuda-toolkit python3-skimage python3-plotly ffmpeg python3-msgpack python3-pandas
>pip install plotly==3.8.1 folium moviepy cupy --upgrade --user

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

**Git over http (Optional)**
On the remote server:
>sudo ln -s /home/exequiel/repositories /var/www/html/
>sudo vim /etc/apache2/conf-available/cgi-enabled.conf
    <Directory "/usr/lib/git-core/">
      Options +ExecCGI
      Require all granted
    </Directory>
    SetEnv GIT_PROJECT_ROOT /var/www/html/repositories
    SetEnv GIT_HTTP_EXPORT_ALL
    ScriptAlias /repositories /usr/lib/git-core/git-http-backend
>sudo service apache2 restart
source:https://www.webmasterworld.com/apache/4711334.htm

**To clone the project**
>git clone exequiel@pluton.hopto.org:/home/exequiel/repositories/aedes_aegypti.git
or
>git clone ssh://exequiel@pluton.hopto.org:8080/home/exequiel/repositories/aedes_aegypti.git
or
>git clone http://pluton/repositories/aedes_aegypti.git

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

**Installing pybind11**
>sudo apt install python3-dev
>git clone https://github.com/pybind/pybind11.git
>sudo cp pybind11/include/pybind11 /usr/include/  -r

**Installing Eigen**
>git clone https://gitlab.com/libeigen/eigen.git
>sudo cp eigen/Eigen /usr/include/  -r

**c++ binding**
>g++ -std=c++17 -Wall -O3 -march=native -shared -fPIC -I/usr/include/python3.6m src/cpp/otero_precipitation_wrapper.cpp -o src/otero_precipitation_wrapper.so

**c++ profiling**
>g++ -std=c++17 -O3 -march=native src/cpp/main.cpp; sudo perf record -g ./a.out;sudo perf report
>g++ -std=c++17 -O3 -march=native src/cpp/main.cpp -pg;  time ./a.out ;gprof -l > gprof2.txt

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

#Strange error
error:
  python: /usr/include/Eigen/src/Core/Block.h:147: Eigen::Block<XprType, BlockRows, BlockCols, InnerPanel>::Block(XprType&, Eigen::Index, Eigen::Index, Eigen::Index, Eigen::Index) [with XprType = Eigen::Array<double, -1, -1>; int BlockRows = 1; int BlockCols = -1; bool InnerPanel = false; Eigen::Index = long int]: Assertion 'startRow >= 0 && blockRows >= 0 && startRow <= xpr.rows() - blockRows && startCol >= 0 && blockCols >= 0 && startCol <= xpr.cols() - blockCols' failed.
  Aborted (core dumped)
Aparently when you import matplotlib, it changes the behaviour of std::stod (wtf!) , so in the pybind module, the configuration::getScalar returned vBS_lh=0 (instead 0.1), so m=0 (instead of 100) and that ended up messing the line mE0(0,Eigen::all)= initial_condition(0)*this->parameters.vBS_d; in otero_precipitation.h
Solution: remove matplotlib
Source: https://github.com/pybind/pybind11/issues/2042
