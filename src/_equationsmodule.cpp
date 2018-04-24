//>g++  -Wall -std=c++11 -I/usr/include/python2.7 -fpic  _equationsmodule.cpp -shared -lboost_python -o _equationsmodule.so -O3
#include "_equationsmodule.h"
#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(_equations)
{
        def("vR_D", vR_D);
        def("dvE", dvE);
        def("dvL", dvL);
        def("dvP", dvP);
        def("dA1", dA1);
        def("dA2", dA2);
        def("dW", dW);

        def("initialize", initialize);

}
