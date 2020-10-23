**Dependencies**
>sudo apt install python3-scipy python3-grib python3-netcdf4 python3-matplotlib python3-gdal python3-line-profiler nvidia-cuda-toolkit python3-skimage python3-plotly ffmpeg python3-msgpack python3-pandas
>pip install plotly==3.8.1 folium moviepy cupy --upgrade --user

**Run**
>python src/tests.py 1



**c++ binding**
>g++ -std=c++17 -Wall -O3 -march=native -shared -fPIC -I/usr/include/python3.6m src/cpp/otero_precipitation_wrapper.cpp -o src/otero_precipitation_wrapper.so

**c++ profiling**
>g++ -std=c++17 -O3 -march=native src/cpp/main.cpp; sudo perf record -g ./a.out;sudo perf report
>g++ -std=c++17 -O3 -march=native src/cpp/main.cpp -pg;  time ./a.out ;gprof -l > gprof2.txt
