#ifndef EquationsH
#define EquationsH
#include <iostream>
#include <math.h>

//TODO:change the name of the following arrays
const double R_D_298K_[]={0.24,0.2088,0.384,0.216,0.372};
const double deltaH_A_[]={10798.0,26018.0,14931.0,15725.0,15725.0};
const double deltaH_H_[]={100000.0,55990.0,-472379.00,1756481.0,1756481.0};
const double T_1_2H_[]={14184.0,304.6,148.0,447.2,447.2};

//example with a class
class Equations {
    static Equations* equations;
    Equations(){}
    public:
    //for Singletons we are using static local, so the destructor is called when the program exits
    static Equations& getInstance(){
        static Equations equations;
        return equations;
    }

    double sum(double a,double b){
        return a+b;
    }

    double mul(double a,double b){
        return a*b;
    }

    double R_D(int stage,double T_t){//day^-1
        double R=1.987; // universal gas constant
        double R_D_298K=R_D_298K_[stage];
        double deltaH_A=deltaH_A_[stage];
        double deltaH_H=deltaH_H_[stage];
        double T_1_2H=T_1_2H_[stage]; // K
        return R_D_298K * (T_t/298.0) * exp( (deltaH_A/R)* ((1.0/298.0)- (1.0/T_t)) ) / ( 1.0+ exp( (deltaH_H/R)* ((1.0/T_1_2H)-(1.0/T_t)) ) );
    }
};

#endif

//example without a class
double sum(double a,double b){
    return a+b;
}

double mul(double a,double b){
  return a*b;
}

double R_D(int stage,double T_t){//day^-1
    double R=1.987; // universal gas constant
    double R_D_298K=R_D_298K_[stage];
    double deltaH_A=deltaH_A_[stage];
    double deltaH_H=deltaH_H_[stage];
    double T_1_2H=T_1_2H_[stage]; // K
    return R_D_298K * (T_t/298.0) * exp( (deltaH_A/R)* ((1.0/298.0)- (1.0/T_t)) ) / ( 1.0+ exp( (deltaH_H/R)* ((1.0/T_1_2H)-(1.0/T_t)) ) );
}

//this main method has no effect on python, is just to test from a c++ environment
int main(){
  std::cout << Equations::getInstance().sum(3,2);
}
