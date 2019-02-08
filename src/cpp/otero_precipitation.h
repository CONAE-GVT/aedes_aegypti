#ifndef ModelH
#define ModelH

#include "types.h"
#include "configuration.h"
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
    std::vector<scalar> time_range;

    Model():Model(Configuration("resources/otero_precipitation.cfg")){}

    Model(Configuration configuration){
        this->parameters.BS_a=configuration.getScalar("breeding_site","amount");
        this->parameters.vBS_h=configuration.getTensor("breeding_site","height");//#in cm
        this->parameters.vBS_s=configuration.getTensor("breeding_site","surface");//#in cm^2
        this->parameters.vBS_d=configuration.getTensor("breeding_site","distribution");//#distribution of BS. Sum must be equals to 1
        this->parameters.vBS_W0=configuration.getTensor("breeding_site","initial_water");
        this->parameters.vBS_mf=configuration.getTensor("breeding_site","manually_filled");//#in percentage of capacity
        this->parameters.n=this->parameters.vBS_d.size();


        unsigned int n=this->parameters.n;
        this->parameters.EGG=std::slice(0,n,1);//#in R^n
        this->parameters.LARVAE=std::slice(n,n,1);//#in R^n
        this->parameters.PUPAE=std::slice(2*n,n,1);//#in R^n
        this->parameters.ADULT1=3*n;//#in R
        this->parameters.ADULT2=3*n+1;//#in R

        this->parameters.vAlpha0=configuration.getTensor("biology","alpha0");//#constant to be fitted

        this->parameters.location=configuration.get("location","name");
        this->start_date=configuration.get("simulation","start_date");
        this->end_date=configuration.get("simulation","end_date");
        this->parameters.initial_condition=configuration.getTensor("simulation","initial_condition");
        std::string WEATHER_DATA_FILENAME="data/public/"+this->parameters.location+".csv";
        this->parameters.weather=Weather(WEATHER_DATA_FILENAME, this->start_date ,this->end_date );
        unsigned int days=Utils::getDaysFromCsv(WEATHER_DATA_FILENAME, this->start_date ,this->end_date );
        for(unsigned int i=0;i<days;i++) this->time_range.push_back(i);


        this->parameters.mf=[](scalar t) { return (1.-std::min(int(t)%7,1))* (sin(2.*M_PI*t + 3.*M_PI/2.) +1.); };//<---- this is implemented different in python
        std::vector<tensor> W = RK::solve(waterEquations,this->parameters.vBS_W0,time_range,this->parameters,20);
        std::vector<std::function<scalar(scalar)>> vWaux=std::vector<std::function<scalar(scalar)>>();
        for(unsigned int j=0;j<n;j++) vWaux.push_back( Utils::getSpline(this->time_range,Utils::getColumn(W,j)) );
        this->parameters.vW=[vWaux](scalar t){tensor vW_t=tensor(vWaux.size()); for(unsigned int i=0;i<vWaux.size();i++) vW_t[i]=std::max(vWaux[i](t),0.0); return vW_t;};

    }


    std::vector<tensor> solveEquations(){
        std::vector<scalar> time_range= this->time_range;
        tensor Y0=this->parameters.initial_condition;
        std::vector<tensor> Y=RK::solve(diff_eqs,Y0,time_range,this->parameters,20);
        return Y;
    }

};

#endif
