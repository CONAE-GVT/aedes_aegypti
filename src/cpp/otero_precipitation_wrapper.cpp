//g++ -std=c++17 -Wall -O3 -march=native -shared -fPIC -I/usr/include/python3.6m src/cpp/otero_precipitation_wrapper.cpp -o src/otero_precipitation_wrapper.so
#include "otero_precipitation.h"
#include "configuration.h"
//order is important! if I put the pybind includes before mine, the result contains nans or infs.(The problem seems to be just with the eigen import)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

//https://pybind11.readthedocs.io/en/stable/classes.html
PYBIND11_MODULE(otero_precipitation_wrapper, m) {
    pybind11::class_<Model>(m, "Model")
        .def(pybind11::init<const std::string>())
        .def("solveEquations", &Model::solveEquations)
        .def_readonly("start_date", &Model::start_date)
        .def_readonly("end_date", &Model::end_date)
        .def_readonly("time_range", &Model::time_range)
        .def_readonly("Y", &Model::Y);
}
