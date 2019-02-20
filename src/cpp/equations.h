#ifndef EquationsH
#define EquationsH

//TODO:make this a class?
#include "types.h"
#include "parameters.h"
#include "utils.h"



tensor vR_D_298K={0.24,0.2088,0.384,0.216,0.372};
//#ro_25_=[0.01066,0.00873,0.01610,0.00898] #replaced by R_D_298K. which:  R_D_298K ~ ro_25*24  #(24 hours)
tensor vDeltaH_A={10798.0,26018.0,14931.0,15725.0,15725.0};
tensor vDeltaH_H={100000.0,55990.0,-472379.00,1756481.0,1756481.0};// #-472379 vs. -473379
tensor vT_1_2H={14184.0,304.6,148.0,447.2,447.2};

//<precipitation related functionality>
//#Ivanov
scalar QR(scalar RH_t,scalar T_t){//#in cm/day
    return 6e-5* std::pow(25 + T_t-273.15, 2) * (100.-RH_t) * 0.1;//#cm
}

scalar QG(scalar p_t){//#Quantity gathered#in cm
    return p_t*0.1;//#cm
}

/*
    { QG(BS_s,t)-QR(BS_s,T(t))    if 0 < W < BS_h
dW= { QG(BS_s,t)                  if W <= 0.0
    { -QR(BS_s,T(t))               if W >= BS_h
Note: in the implementation we needed to add functions to make function continuous, otherwise odeint breaks
*/
scalar dW(scalar W, scalar BS_h, scalar T_t, scalar p_t, scalar RH_t){//#in cm/day
    scalar epsilon=1e-1;//#1mm
    if(0+epsilon < W && W < BS_h-epsilon)
        return QG(p_t)-QR(RH_t,T_t);
    else if(W <= 0.0+epsilon)
        return QG(p_t) - QR(RH_t,T_t)*(W/epsilon);
    else// if( W >= BS_h-epsilon)//commented out to avoid warning(should have no effect)
        return QG(p_t)*((BS_h-W)/epsilon) - QR(RH_t,T_t);
}

scalar a0(scalar W){
    return 70.0* W;
}

scalar gamma(scalar L, scalar BS, scalar W);
tensor vGamma(const tensor& vL, tensor vBS_a, const tensor& vW){
    tensor vGamma_t=tensor(vW.size());
    for(unsigned int i=0;i<vL.size();i++) vGamma_t[i]=gamma(vL[i],vBS_a[i],vW[i]);
    return vGamma_t;
}

tensor waterEquations(const tensor& vW,scalar t, Parameters& parameters){
    scalar T_t=parameters.weather.T(t);
    scalar p_t=parameters.weather.p(t);
    scalar RH_t=parameters.weather.RH(t);
    tensor vmf_t=parameters.mf(t)*parameters.vBS_mf*parameters.vBS_h*10.;//#% -> cm -> mm
    tensor vBS_h=parameters.vBS_h;
    tensor dWdt=tensor(vW.size());
    for(unsigned int i=0;i<vW.size();i++) dWdt[i]=dW(vW[i],vBS_h[i],T_t,p_t+vmf_t[i],RH_t);
    return dWdt;
}
//</precipitation related functionality v>

tensor vR_D(scalar T_t){//#day^-1
    scalar R=1.987; //# universal gas constant
    return vR_D_298K * (T_t/298.0) * std::exp( (vDeltaH_A/R)* ((1.0/298.0)- (1.0/T_t)) ) / ( 1.0+ std::exp( (vDeltaH_H/R)* ((1.0/vT_1_2H)-(1.0/T_t)) ) );
}

scalar gamma(scalar L, scalar BS, scalar W){
    scalar epsilon=1e-4;
    if(BS==0 || W <epsilon)//#W *1000./BS_s <0.1
        return 1.0;//#no water total inhibition#1960 Aedes aegypti (L.), The Yellow Fever Mosquito(Page 165)
    if(L/BS<=a0(W)-epsilon)
        return  0;
    else if(a0(W)-epsilon < L/BS && L/BS<a0(W)+epsilon){
        //#a (a0-e) + b=0 => b=-a (a0 -e)
        //#a (a0 + e) + b=0.63 => a(a0+e) - a(a0-e) = 2 a e = 0.63 =>a=0.63/(2 e)
        scalar a=0.63/(2.0*epsilon);
        scalar b=-a*(a0(W)-epsilon);
        return a * (L/BS) + b;
    }
    else// if(L/BS>=a0(W)+epsilon)//commented out to avoid warning(should have no effect)
        return 0.63;
}

matrix ovsp(const tensor& vW,const tensor& vBS_d,const tensor& vW_l,const matrix& mBS_l){
    scalar epsilon=1e-4;
    tensor vf=vW/(vW+epsilon) * vBS_d;
    if(vf.max()>epsilon) vf=vf/vf.sum();
    matrix ovsp_t=matrix(tensor(vW_l.size()),mBS_l.size());
    for(unsigned int i=0;i<mBS_l.size();i++)
        for(unsigned int j=0;j<mBS_l[i].size();j++)
            if(mBS_l[i][j]==std::floor(vW_l[j])) ovsp_t[i][j]=1.;
            else ovsp_t[i][j]=0.;
    return ovsp_t*vf;//#TODO: check this
}

matrix wetMask(const tensor& vW_l, const matrix& mBS_l){
    matrix mask=matrix(tensor(vW_l.size()),mBS_l.size());
    for(unsigned int i=0;i<mBS_l.size();i++)
        for(unsigned int j=0;j<mBS_l[i].size();j++)
            if(mBS_l[i][j]<=vW_l[j]) mask[i][j]=1.;
            else mask[i][j]=0.;
    return mask;
}

matrix dvE(const matrix& mE,const tensor& vL,scalar A1,scalar A2,const tensor& vW_t, scalar BS_a,const tensor&  vBS_d,scalar elr,scalar ovr1, scalar ovr2,const matrix& wet_mask,const tensor& vW_l,const matrix& mBS_l){
    scalar egn=63.0;
    scalar me=0.01;//#mortality of the egg, for T in [278,303]
    matrix ovsp_t=ovsp(vW_t,vBS_d,vW_l,mBS_l);
    return egn  *( ovr1 *A1  + ovr2* A2)*ovsp_t - me * mE - tensor(elr* (1.-vGamma(vL,BS_a*vBS_d,vW_t))) * mE*wet_mask;
}

tensor dvL(const matrix& mE,const tensor& vL,const tensor& vW,scalar T_t,scalar BS_a,const tensor& vBS_d,scalar elr,scalar lpr,const tensor& vAlpha0,const matrix& wet_mask){
    scalar ml=0.01 + 0.9725 * std::exp(-(T_t-278.0)/2.7035);//#mortality of the larvae, for T in [278,303]
    scalar mdl=2.;//#mortality of dry larvae.TODO:Unjustified!
    tensor vAlpha=vAlpha0/(BS_a*vBS_d);
    scalar epsilon=1e-4;
    return elr* (1.-vGamma(vL,BS_a*vBS_d,vW)) * Utils::sumAxis0(mE*wet_mask) - ml*vL - vAlpha* vL*vL - lpr *vL - mdl*(1.- vW/(vW+epsilon))*vL;
}

tensor dvP(const tensor& vL,const tensor& vP,scalar T_t, scalar lpr,scalar par){
    scalar mp=0.01 + 0.9725 * std::exp(-(T_t-278.0)/2.7035);//#death of pupae
    return lpr*vL - mp*vP  - par*vP;
}

scalar dA1(const tensor& vP,scalar A1, scalar par,scalar ovr1){
    scalar ef=0.83;//#emergence factor
    scalar ma=0.09;//#for T in [278,303]
    return (par*ef*vP/2.0).sum() - ma*A1 - ovr1*A1;
}

scalar dA2(scalar A1,scalar A2,scalar ovr1){
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
    scalar BS_lh=parameters.BS_lh;
    tensor vBS_d=parameters.vBS_d;
    tensor vAlpha0=parameters.vAlpha0;
    unsigned int m=parameters.m;
    unsigned int n=parameters.n;
    matrix mBS_l=parameters.mBS_l;

    tensor vW_t=parameters.vW(t);
    matrix vE=Utils::tensorToMatrix(Y[parameters.EGG],m,n);//equivalent to reshape and transpose
    tensor vL=Y[parameters.LARVAE];
    tensor vP=Y[parameters.PUPAE];
    scalar A1=Y[parameters.ADULT1];
    scalar A2=Y[parameters.ADULT2];
    tensor vW_l=vW_t/BS_lh;
    matrix wet_mask=wetMask(vW_l,mBS_l);

    tensor dY=tensor(Y.size());
    dY[parameters.EGG]    = Utils::matrixToTensor(dvE(vE,vL,A1,A2,vW_t,BS_a,vBS_d,elr,ovr1,ovr2,wet_mask,vW_l,mBS_l));
    dY[parameters.LARVAE] = dvL(vE,vL,vW_t,T_t,      BS_a,vBS_d,elr,lpr,vAlpha0,wet_mask);
    dY[parameters.PUPAE]  = dvP(vL,vP,T_t,lpr,par);
    dY[parameters.ADULT1] = dA1(vP,A1,par,ovr1);
    dY[parameters.ADULT2] = dA2(A1,A2,ovr1);

    return dY;//   # For odeint
}

#endif
