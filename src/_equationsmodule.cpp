//>g++  -Wall -std=c++11 -I/usr/include/python2.7 -fpic  _equationsmodule.cpp -shared -lboost_python -o _equationsmodule.so -O3
#include "_equationsmodule.h"
#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(_equations)
{
    def("getInstance", Equations::getInstance,
        return_value_policy<reference_existing_object>());
    class_<Equations>("Equations",no_init)
        .def("sum",&Equations::sum)
        .def("mul",&Equations::mul)
        .def("R_D",&Equations::R_D);

        def("sum", sum);
        def("mul", mul);
        def("R_D", R_D);
}
