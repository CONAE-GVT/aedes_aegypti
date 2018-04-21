#ifndef EquationsH
#define EquationsH
#include <iostream>
#include <math.h>
#include <vector>

//TODO:change the name of the following arrays
const double R_D_298K_[]={0.24,0.2088,0.384,0.216,0.372};
const double deltaH_A_[]={10798.0,26018.0,14931.0,15725.0,15725.0};
const double deltaH_H_[]={100000.0,55990.0,-472379.00,1756481.0,1756481.0};
const double T_1_2H_[]={14184.0,304.6,148.0,447.2,447.2};

const int EGG=0;
const int LARVAE=1;
const int PUPAE=2;
const int ADULT1=3;
const int ADULT2=4;



double np_dot(std::vector<double> v, std::vector<double> w){//TODO: use std::inner_product(std::begin(a), std::end(a), std::begin(b), 0.0);
    if(v.size()!=w.size()) std::cout <<"np_dot Warning: dimensions don't match!!!";
    double a=0;
    for(uint i=0;i<v.size();i++){
        a+=v[i]*w[i];
    }
    return a;
}
std::vector<double> np_sum(std::vector<double> v, std::vector<double> w){//TODO:use algo
    if(v.size()!=w.size()) std::cout <<"np_sum Warning: dimensions don't match!!!";
    for(uint i=0;i<v.size();i++){
        v[i]+=w[i];
    }
    return v;
}
std::vector<double> np_prod(double a, std::vector<double> v){//TODO:use algo
    for(uint i=0;i<v.size();i++){
        v[i]*=a;
    }
    return v;
}

std::vector<double> np_ones(int n){//TODO:use algo
    return std::vector<double>(n,1.0);
}

/************************************************/
//<precipitation related functionality v>


//NOTE:EPA evaporation in Technical guidance for hazards analysis, eq (7) page G-3
//TODO:check QS for spilled water
double QR(double u,double BS_s,double T_t){//in l/day
    double U=u * 1000.0/(3600.0);//#step(wind_speed,t) * 1000.0/(60.0*60.0)#in m/s #km/h->m/s
    double MW=18.01528;//#molecular wight of water in g/mol (gmol often called mol src:https://chemistry.stackexchange.com/questions/53455/g-gmol-vs-lb-lbmol)
    double A=BS_s * 0.00107;//#in ft^2 #cm^2->ft^2
    double VP=pow(10.0,(8.07131-(1730.63/(233.426+T_t-273.15))));//#Vapor pressure by Antoine Equation in mmHg
    double R=82.05; //#in atm cm 3 /g mole
    return ( (0.106 * pow(U,0.78) * pow(MW,(2.0/3.0))* A * VP)/(R* T_t) ) * 453.59237/(1.0/(60.0*24.0)) *1./1000.; //#Rate of release to air ml/day #(Ibs/min) ->  ml/day ->l/day
}

double  QG(double BS_s,double p_t,double t){//#Quantity gathered#in litres
    return (BS_s * p_t*0.1) * 1.0 * 1./1000.;//#*1cm^3=1ml -> l
}

/*
    { QG(BS_s,t)-QR(BS_s,T(t))    if 0 < W < BS_c
dW= { QG(BS_s,t)                  if W <= 0.0
    { -QR(BS_s,T(t))               if W >= BS_c
Note: in the implementation we needed to add functions to make function continuous, otherwise odeint breaks
*/

double  dW(double W,double BS_c,double BS_s,double T_t,double p_t,double wss_t,double t){//#in l/day
    double epsilon=1e-3;
    if(0+epsilon < W && W < BS_c-epsilon){
        return QG(BS_s,p_t,t)-QR(wss_t,BS_s,T_t);
    }else if(W <= 0.0+epsilon){
        return QG(BS_s,p_t,t) - QR(wss_t,BS_s,T_t)*(W/epsilon);
    }else if( W >= BS_c-epsilon){
        return QG(BS_s,p_t,t)*((BS_c-W)/epsilon) - QR(wss_t,BS_s,T_t);
    }
    return 0;//to avoid a warning
}
double a0(double W){
    return 70.0* W;
}

//#</precipitation related functionality v>

double R_D(int stage,double T_t){//day^-1
    double R=1.987; // universal gas constant
    double R_D_298K=R_D_298K_[stage];
    double deltaH_A=deltaH_A_[stage];
    double deltaH_H=deltaH_H_[stage];
    double T_1_2H=T_1_2H_[stage]; // K
    return R_D_298K * (T_t/298.0) * exp( (deltaH_A/R)* ((1.0/298.0)- (1.0/T_t)) ) / ( 1.0+ exp( (deltaH_H/R)* ((1.0/T_1_2H)-(1.0/T_t)) ) );
}


double  gamma(double L,double BS,double W){
    double epsilon=1e-4;
    if(BS==0 || W <epsilon){//#W *1000./BS_s <0.1
        return 1.0;//#no water total inhibition#1960 Aedes aegypti (L.), The Yellow Fever Mosquito(Page 165)
    }
    if(L/BS<=a0(W)-epsilon){
        return 0;
    }
    else if(a0(W)-epsilon < L/BS && L/BS <a0(W)+epsilon){
        //#a (a0-e) + b=0 => b=-a (a0 -e)
        //#a (a0 + e) + b=0.63 => a(a0+e) - a(a0-e) = 2 a e = 0.63 =>a=0.63/(2 e)
        double a=0.63/(2.0*epsilon);
        double b=-a*(a0(W)-epsilon);
        return a * (L/BS) + b;
    }
    else if(L/BS>=a0(W)+epsilon){
        return 0.63;
    }
    return 0;//to avoid a warning
}

std::vector<double> vGamma(std::vector<double> vL,std::vector<double> vBS_a,std::vector<double> vW){
    std::vector<double> result=std::vector<double>();
    for(uint i=0;i<vL.size();i++){
        result.push_back(gamma(vL[i],vBS_a[i],vW[i]) );
    }
    return result;
}



double beta(std::vector<double> vW,std::vector<double> vBS_od,std::vector<double> vBS_id){//#~1 if we have water,~0 if we dont
    double result=0;
    for(uint i=0; i<vBS_od.size();i++){
        result+= vBS_od[i]*vW[i]/(vW[i]+1e-4);
    }
    for(uint i=0; i<vBS_id.size();i++){
     result+=vBS_id[i];//#TODO:check!
    }
    return result;
}

double dE(double E,double L,double A1,double A2,std::vector<double> vW,double T_t,double BS_a,std::vector<double> vBS_ic,std::vector<double> vBS_od,std::vector<double> vBS_id,int n,int m){
    double egn=63.0*beta(vW,vBS_od,vBS_id);//#The amount of eggs goes to zero when vW goes to zero.In the 1-dimensional case. when W=1e-3, egn(1e-3)~egn(0.5)/2#TODO:instead of 1e-3, it whould be a value related to the min water needed to lay eggs
    double me=0.01;//#mortality of the egg, for T in [278,303]
    double elr=R_D(EGG,T_t);
    double ovr1=R_D(ADULT1,T_t);
    double ovr2=R_D(ADULT2,T_t);
    std::vector<double> v1=np_ones(n);
    //#  ( (1,1,..1) - vGamma(vL,vBS_a,vW) ) . vE = (1- γ(vL[1],vBS_a[1],vW[1],...) ) . vE= Σ (1- γ(vL[i],vBS_a[i],vW[i])) * vE[i]
    std::vector<double> v1_vGamma=np_sum(v1, np_prod(-1.,vGamma(np_prod(L,vBS_od),np_prod(BS_a,vBS_od),vW    )) );
    double inh_o=np_dot( v1_vGamma, np_prod(E,vBS_od) );
    v1=np_ones(m);
    //#  ( (1,1,..1) - vGamma(vL,vBS_a,vBS_ic) ) . vE = (1-γ(vL[1],vBS_a[1],vBS_ic[1]) ,... ) . vE= Σ (1- γ(vL[i],vBS_a[i],vBS_ic[i])) * vE[i]
    v1_vGamma=np_sum(v1, np_prod(-1.,vGamma(np_prod(L,vBS_id),np_prod(BS_a,vBS_id),vBS_ic    )) );
    double inh_i=np_dot(v1_vGamma ,np_prod(E,vBS_id) );
    return egn*( ovr1 *A1  + ovr2* A2) - me *E - elr* (inh_o + inh_i );
}
double dL(double E,double L,std::vector<double> vW,double T_t,double BS_a,std::vector<double> vBS_ic,std::vector<double> vBS_od,std::vector<double> vBS_id,int n,int m){
    double elr=R_D(EGG,T_t);
    double lpr=R_D(LARVAE,T_t);
    double ml=0.01 + 0.9725 * exp(-(T_t-278.0)/2.7035);//#mortality of the larvae, for T in [278,303];
    double alpha0=1.5;//#Parameter to be fitted #1.0#HARDCODED!!!
    double alpha=alpha0/BS_a;//#Σ vα[i] * vL[i]^2= Σ α0/(BS_a* vBS_d[i]) * (L*vBS_d[i])^2 = Σ α0/BS_a * L^2 * vBS_d[i] = α *L^2 *Σ BS_d[i]=α L^2 #Note: on why this is ok even though BS changed.
    std::vector<double> v1=np_ones(n);
    std::vector<double> v1_vGamma=np_sum(v1, np_prod(-1.,vGamma(np_prod(L,vBS_od),np_prod(BS_a,vBS_od),vW    )) );
    double inh_o=np_dot( v1_vGamma, np_prod(E,vBS_od) );
    v1=np_ones(m);
    v1_vGamma=np_sum(v1, np_prod(-1.,vGamma(np_prod(L,vBS_id),np_prod(BS_a,vBS_id),vBS_ic    )) );
    double inh_i=np_dot(v1_vGamma ,np_prod(E,vBS_id) );
    return elr* (inh_o+inh_i ) - ml*L - alpha* L*L - lpr *L -35.6464*(1.-beta(vW,vBS_od,vBS_id))*L;//#-24.*(1.-beta(vW))*L# -log(1e-4/5502.)/(1.)=17.823207313460703
}
double dP(double L,double P,double T_t){
    double lpr=R_D(LARVAE,T_t);
    double par=R_D(PUPAE,T_t);
    double mp=0.01 + 0.9725 * exp(-(T_t-278.0)/2.7035);//#death of pupae
    return lpr*L - mp*P  - par*P;
}

double dA1(double P,double A1,double T_t){
    double par=R_D(PUPAE,T_t);
    double ovr1=R_D(ADULT1,T_t);
    double ef=0.83;//#emergence factor
    double ma=0.09;//#for T in [278,303]
    return par*ef*P/2.0 - ma*A1 - ovr1*A1;
}

double dA2(double A1,double A2,double T_t){
    double ovr1=R_D(ADULT1,T_t);
    double ma=0.09;//#for T in [278,303]
    return ovr1*A1 - ma*A2;
}

/************************************************/
//this main method has no effect on python, is just to test from a c++ environment


#endif
int main(){
  std::cout << dA2(1,2,3);
}
