#ifndef ModelH
#define ModelH

#include "parameters.h"
#include "weather.h"
#include "equations.h"
#include "rk.h"

class Model
{
    public:
    Parameters parameters;
    std::string start_date;
    std::string end_date;
    std::vector<double> time_range;

    Model(){
        this->parameters.BS_a=50.;
        this->parameters.vBS_h={10.,10.};//#in cm
        this->parameters.vBS_s={50,50};//#in cm^2
        this->parameters.vBS_d={0.5,0.5};//#distribution of BS. Sum must be equals to 1
        this->parameters.vBS_W0={0,0};
        this->parameters.vBS_mf={0.,0.};//#in percentage of capacity
        this->parameters.n=this->parameters.vBS_d.size();

        //TODO:implement
        unsigned int n=this->parameters.n;
        this->parameters.EGG=std::slice(0,n,1);//#in R^n
        this->parameters.LARVAE=std::slice(n,n,1);//#in R^n
        this->parameters.PUPAE=std::slice(2*n,n,1);//#in R^n
        this->parameters.ADULT1=3*n;//#in R
        this->parameters.ADULT2=3*n+1;//#in R

        this->parameters.vAlpha0={1.5,1.5};//#constant to be fitted

        this->start_date="2015-7-1";
        this->end_date="2019-1-15";
        this->parameters.initial_condition={100.,100., 0,0, 0,0, 0, 0};
        this->parameters.weather=Weather("data/public/wunderground.csv", this->start_date ,this->end_date );
        double h=1/2.;
        unsigned int days=Utils::getDaysFromCsv("data/public/wunderground.csv", this->start_date ,this->end_date );
        for(unsigned int i=0;i<days/h;i++) this->time_range.push_back(i*h);

        //TODO:implement water!!!
    }


    std::vector<std::valarray<double>> solveEquations(){
        std::vector<double> time_range= this->time_range;
        std::valarray<double> Y0=this->parameters.initial_condition;
        std::vector<std::valarray<double>> Y=RK::solve(diff_eqs,Y0,time_range,this->parameters,20);
        return Y;
    }

};

#endif
