#ifndef EquationsH
#define EquationsH

//TODO:make this a class?
#include "types.h"
#include "parameters.h"



tensor vR_D_298K={0.24,0.2088,0.384,0.216,0.372};
//#ro_25_=[0.01066,0.00873,0.01610,0.00898] #replaced by R_D_298K. which:  R_D_298K ~ ro_25*24  #(24 hours)
tensor vDeltaH_A={10798.0,26018.0,14931.0,15725.0,15725.0};
tensor vDeltaH_H={100000.0,55990.0,-472379.00,1756481.0,1756481.0};// #-472379 vs. -473379
tensor vT_1_2H={14184.0,304.6,148.0,447.2,447.2};

tensor vR_D(scalar T_t){//#day^-1
    scalar R=1.987; //# universal gas constant
    return vR_D_298K * (T_t/298.0) * std::exp( (vDeltaH_A/R)* ((1.0/298.0)- (1.0/T_t)) ) / ( 1.0+ std::exp( (vDeltaH_H/R)* ((1.0/vT_1_2H)-(1.0/T_t)) ) );
}

//TODO:IMPLEMENT!!!
tensor vGamma(const tensor& vL,tensor vBS_a,const tensor& vW){
    return tensor(vW.size());
}
tensor ovsp(const tensor& vBS_d,const tensor& vW){
    scalar epsilon=1e-4;
    tensor vf=vW/(vW+epsilon) * vBS_d;
    if(vf.max()<1e-20) return vf;
    else return vf/vf.sum();//#TODO: check this
}

tensor dvE(const tensor& vE,const tensor& vL,scalar A1,scalar A2,const tensor& vW, scalar T_t, scalar BS_a,const tensor&  vBS_d,scalar elr,scalar ovr1, scalar ovr2){
    scalar egn=63.0;
    scalar me=0.01;//#mortality of the egg, for T in [278,303]
    return egn*( ovr1 *A1  + ovr2* A2)*ovsp(vW,vBS_d) - me * vE - elr* (1.-vGamma(vL,BS_a*vBS_d,vW)) * vE;
}

tensor dvL(const tensor& vE,const tensor& vL,const tensor& vW,scalar T_t,scalar BS_a,const tensor& vBS_d,scalar elr,scalar lpr,const tensor& vAlpha0){
    scalar ml=0.01 + 0.9725 * std::exp(-(T_t-278.0)/2.7035);//#mortality of the larvae, for T in [278,303]
    tensor vAlpha=vAlpha0/(BS_a*vBS_d);
    return elr* (1.-vGamma(vL,BS_a*vBS_d,vW)) * vE - ml*vL - vAlpha* vL*vL - lpr *vL ;//#-35.6464*(1.-beta(vW,vBS_od,vBS_id))*L#-24.*(1.-beta(vW))*L# -log(1e-4/5502.)/(1.)=17.823207313460703
}

tensor dvP(const tensor& vL,const tensor& vP,scalar T_t, scalar lpr,scalar par){
    scalar mp=0.01 + 0.9725 * std::exp(-(T_t-278.0)/2.7035);//#death of pupae
    return lpr*vL - mp*vP  - par*vP;
}

scalar dA1(const tensor& vP,scalar A1,scalar T_t, scalar par,scalar ovr1){
    scalar ef=0.83;//#emergence factor
    scalar ma=0.09;//#for T in [278,303]
    return (par*ef*vP/2.0).sum() - ma*A1 - ovr1*A1;
}

scalar dA2(scalar A1,scalar A2,scalar T_t,scalar ovr1){
    scalar ma=0.09;//#for T in [278,303]
    return ovr1*A1 - ma*A2;
}


tensor diff_eqs(const tensor& Y,scalar t,Parameters& parameters){
    scalar T_t=parameters.weather.T(t);
    tensor vR_D_t=vR_D(T_t);
    scalar elr=vR_D_t[0];
    scalar lpr=vR_D_t[1];
    scalar par=vR_D_t[2];
    scalar ovr1=vR_D_t[3];
    scalar ovr2=vR_D_t[4];

    scalar BS_a=parameters.BS_a;
    tensor vBS_d=parameters.vBS_d;
    tensor vAlpha0=parameters.vAlpha0;

    tensor vW_t=tensor(1,parameters.ADULT2);//parameters.vW(t)TODO:implement
    tensor vE=Y[parameters.EGG];
    tensor vL=Y[parameters.LARVAE];
    tensor vP=Y[parameters.PUPAE];
    scalar A1=Y[parameters.ADULT1];
    scalar A2=Y[parameters.ADULT2];

    tensor dY=tensor(Y.size());
    dY[parameters.EGG]    = dvE(vE,vL,A1,A2,vW_t,T_t,BS_a,vBS_d,elr,ovr1,ovr2);
    dY[parameters.LARVAE] = dvL(vE,vL,vW_t,T_t,      BS_a,vBS_d,elr,lpr,vAlpha0);
    dY[parameters.PUPAE]  = dvP(vL,vP,T_t,lpr,par);
    dY[parameters.ADULT1] = dA1(vP,A1,T_t,par,ovr1);
    dY[parameters.ADULT2] = dA2(A1,A2,T_t,ovr1);

    return dY;//   # For odeint
}

#endif
