//g++ -O3 -Wall -I/usr/include/python3.6m/   -fpic  src/cpp/_otero_precipitation.cpp -shared  -lboost_python3 -lpython3.6m -o src/_otero_precipitation.so
#include <boost/python.hpp>
#include "otero_precipitation.h"

BOOST_PYTHON_MODULE(_otero_precipitation)
{
    boost::python::class_<Model>("Model").def("solveEquations", &Model::solveEquations);
}
