//g++ -std=c++17 -O3 -march=native -Wall -I/usr/include/python3.6m/   -fpic  src/cpp/otero_precipitation_wrapper.cpp -shared  -lboost_python-py36  -o src/otero_precipitation_wrapper.so
#include <boost/python.hpp>
#include "otero_precipitation.h"
#include "configuration.h"

class ModelWrapper
{
    public:
    Model model;
    ModelWrapper(std::string filename): model(Configuration(filename)){//we initialize this->model this way to avoid default constructor
    }

    boost::python::list solveEquations(){
        std::vector<tensor> Y =this->model.solveEquations();
        boost::python::list _Y=boost::python::list();
        for(unsigned int i=0;i<Y.size();i++){
            boost::python::list _Y_i=boost::python::list();
            for(unsigned int j=0;j<Y[i].size();j++) _Y_i.append(Y[i][j]);
            _Y.append(_Y_i);
        }
        return _Y;
    }
};

BOOST_PYTHON_MODULE(otero_precipitation_wrapper)
{
    boost::python::class_<ModelWrapper>("ModelWrapper",boost::python::init<std::string>()).def("solveEquations", &ModelWrapper::solveEquations);
}
