#ifndef RKH
#define RKH

#include <vector>
#include <valarray>

typedef std::valarray<double> state_type;

class RK{
  public:
    static std::vector<state_type> solve(state_type (*dYdt)(const state_type&, double),state_type& Y0,std::vector<double>& time_range, int steps){//TODO:use c++ functor instead
        //main
        std::vector<state_type> Y=std::vector<state_type>();
        //Y.reserve(time_range.size());//TODO:check if we need this.
        Y.push_back(Y0);//Y[0]=Y0<---initial conditions

        state_type Y_j=Y0;
        for(unsigned int i=0; i<time_range.size()-1;i++){
            double t=time_range[i];
            double h=time_range[i+1]-time_range[i];
            double h_j=h/double(steps);
            for(int j=0;j<steps;j++){
                //#Runge-Kutta's terms
                state_type K_n1=dYdt(Y_j,t);
                state_type K_n2=dYdt(Y_j + (h_j/2.)*K_n1, t + h_j/2.);
                state_type K_n3=dYdt(Y_j + (h_j/2.)*K_n2, t + h_j/2.);
                state_type K_n4=dYdt(Y_j + h_j*K_n3, t + h_j);

                Y_j = Y_j + (h_j/6.0)*(K_n1 + 2.0*K_n2 + 2.0*K_n3 + K_n4);
                t=t+h_j;
            }
            Y.push_back(Y_j);
        }

        return Y;
    }
};

#endif
