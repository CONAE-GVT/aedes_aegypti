#ifndef EquationsH
#define EquationsH

//TODO:make this a class?
#include <valarray>
#include "parameters.h"

typedef std::valarray<double> state_type;//TODO:duplicated in RK

std::valarray<double> vR_D_298K={0.24,0.2088,0.384,0.216,0.372};
//#ro_25_=[0.01066,0.00873,0.01610,0.00898] #replaced by R_D_298K. which:  R_D_298K ~ ro_25*24  #(24 hours)
std::valarray<double> vDeltaH_A={10798.0,26018.0,14931.0,15725.0,15725.0};
std::valarray<double> vDeltaH_H={100000.0,55990.0,-472379.00,1756481.0,1756481.0};// #-472379 vs. -473379
std::valarray<double> vT_1_2H={14184.0,304.6,148.0,447.2,447.2};

std::valarray<double> vR_D(double T_t){//#day^-1
    double R=1.987; //# universal gas constant
    return vR_D_298K * (T_t/298.0) * std::exp( (vDeltaH_A/R)* ((1.0/298.0)- (1.0/T_t)) ) / ( 1.0+ std::exp( (vDeltaH_H/R)* ((1.0/vT_1_2H)-(1.0/T_t)) ) );
}

//TODO:IMPLEMENT!!!
std::valarray<double> vGamma(std::valarray<double> vL,std::valarray<double> vBS_a,std::valarray<double> vW){
    return std::valarray<double>(vW.size());
}
std::valarray<double> ovsp(std::valarray<double> vBS_d,std::valarray<double> vW){
    double epsilon=1e-4;
    std::valarray<double> vf=vW/(vW+epsilon) * vBS_d;
    if(vf.max()<1e-20) return vf;
    else return vf/vf.sum();//#TODO: check this
}

std::valarray<double> dvE(std::valarray<double> vE,std::valarray<double> vL,double A1,double A2,std::valarray<double> vW, double T_t, double BS_a,std::valarray<double>  vBS_d,double elr,double ovr1, double ovr2){
    double egn=63.0;
    double me=0.01;//#mortality of the egg, for T in [278,303]
    return egn*( ovr1 *A1  + ovr2* A2)*ovsp(vW,vBS_d) - me * vE - elr* (1.-vGamma(vL,BS_a*vBS_d,vW)) * vE;
}

std::valarray<double> dvL(std::valarray<double> vE,std::valarray<double> vL,std::valarray<double> vW,double T_t,double BS_a,std::valarray<double> vBS_d,double elr,double lpr,std::valarray<double> vAlpha0){
    double ml=0.01 + 0.9725 * std::exp(-(T_t-278.0)/2.7035);//#mortality of the larvae, for T in [278,303]
    std::valarray<double> vAlpha=vAlpha0/(BS_a*vBS_d);
    return elr* (1.-vGamma(vL,BS_a*vBS_d,vW)) * vE - ml*vL - vAlpha* vL*vL - lpr *vL ;//#-35.6464*(1.-beta(vW,vBS_od,vBS_id))*L#-24.*(1.-beta(vW))*L# -log(1e-4/5502.)/(1.)=17.823207313460703
}

std::valarray<double> dvP(std::valarray<double> vL,std::valarray<double> vP,double T_t, double lpr,double par){
    double mp=0.01 + 0.9725 * std::exp(-(T_t-278.0)/2.7035);//#death of pupae
    return lpr*vL - mp*vP  - par*vP;
}

double dA1(std::valarray<double> vP,double A1,double T_t, double par,double ovr1){
    double ef=0.83;//#emergence factor
    double ma=0.09;//#for T in [278,303]
    return (par*ef*vP/2.0).sum() - ma*A1 - ovr1*A1;
}

double dA2(double A1,double A2,double T_t,double ovr1){
    double ma=0.09;//#for T in [278,303]
    return ovr1*A1 - ma*A2;
}


state_type diff_eqs(const state_type& Y,double t,Parameters& parameters){
    double T_t=parameters.weather.T(t);
    std::valarray<double> vR_D_t=vR_D(T_t);
    double elr=vR_D_t[0];
    double lpr=vR_D_t[1];
    double par=vR_D_t[2];
    double ovr1=vR_D_t[3];
    double ovr2=vR_D_t[4];

    double BS_a=parameters.BS_a;
    std::valarray<double> vBS_d=parameters.vBS_d;
    std::valarray<double> vAlpha0=parameters.vAlpha0;
    unsigned int n=parameters.n;

    std::valarray<double> vW_t=std::valarray<double>(1,parameters.ADULT2);//parameters.vW(t)TODO:implement
    std::valarray<double> vE=Y[parameters.EGG];
    std::valarray<double> vL=Y[parameters.LARVAE];
    std::valarray<double> vP=Y[parameters.PUPAE];
    double A1=Y[parameters.ADULT1];
    double A2=Y[parameters.ADULT2];

    std::valarray<double> dY=std::valarray<double>(Y.size());
    dY[parameters.EGG]    = dvE(vE,vL,A1,A2,vW_t,T_t,BS_a,vBS_d,elr,ovr1,ovr2);
    dY[parameters.LARVAE] = dvL(vE,vL,vW_t,T_t,      BS_a,vBS_d,elr,lpr,vAlpha0);
    dY[parameters.PUPAE]  = dvP(vL,vP,T_t,lpr,par);
    dY[parameters.ADULT1] = dA1(vP,A1,T_t,par,ovr1);
    dY[parameters.ADULT2] = dA2(A1,A2,T_t,ovr1);

    return dY;//   # For odeint
}

#endif
