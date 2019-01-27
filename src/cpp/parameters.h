#ifndef ParametersH
#define ParametersH
#include <valarray>
#include "weather.h"

struct Parameters{
    //config
    double BS_a;
    std::valarray<double> vBS_d;
    std::valarray<double>vBS_s;
    std::valarray<double>vBS_h;
    std::valarray<double>vBS_W0;
    std::valarray<double>vBS_mf;

    std::string start_date;
    std::string end_date;
    std::valarray<double> initial_condition;

    std::valarray<double> vAlpha0;
    //dynamic
    unsigned int n;
    Weather weather;

};

#endif
