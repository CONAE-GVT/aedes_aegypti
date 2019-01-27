#ifndef EquationsH
#define EquationsH

//TODO:make this a class?
#include <valarray>

typedef std::valarray<double> state_type;//TODO:duplicated in RK

state_type diff_eqs( const state_type& Y , double t ){
    const double sigma = 10.0;
    const double R = 28.0;
    const double b = 8.0 / 3.0;
    state_type dY={
                    sigma * ( Y[1] - Y[0] ),
                    R * Y[0] - Y[1] - Y[0] * Y[2],
                    -b * Y[2] + Y[0] * Y[1]
                };

    return dY;
}

#endif
