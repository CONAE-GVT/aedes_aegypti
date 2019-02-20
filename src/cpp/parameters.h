#ifndef ParametersH
#define ParametersH

#include <functional>
#include "types.h"
#include "weather.h"

struct Parameters{
    //config
    scalar BS_a;
    scalar BS_lh;
    tensor vBS_d;
    tensor vBS_s;
    tensor vBS_h;
    tensor vBS_W0;
    tensor vBS_mf;
    matrix mBS_l;

    std::string location;

    std::string start_date;
    std::string end_date;
    tensor initial_condition;

    tensor vAlpha0;
    //dynamic
    unsigned int m;
    unsigned int n;
    std::slice EGG;
    std::slice LARVAE;
    std::slice PUPAE;
    unsigned int ADULT1;
    unsigned int ADULT2;
    Weather weather;
    std::function<scalar(scalar)> mf;
    std::function<tensor(scalar)> vW;

};

#endif
