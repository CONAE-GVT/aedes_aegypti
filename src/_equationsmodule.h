#ifndef EquationsH
#define EquationsH
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
namespace p = boost::python;
namespace np = boost::python::numpy;

/******/
np::ndarray np_exp(np::ndarray v){
  np::ndarray w = np::zeros(v.get_nd(), v.get_shape(), v.get_dtype());
  np::multi_iter iter = np::make_multi_iter(v,w);
  while (iter.not_done())
  {
    double* argument = reinterpret_cast<double*>(iter.get_data(0));
    double* result = reinterpret_cast<double*>(iter.get_data(1));
    *result=exp(*argument);
    iter.next();
  }
  return w;
}
np::ndarray np_exp(boost::python::api::object v){
  return np_exp(np::from_object(v));
}


double np_dot(np::ndarray v,np::ndarray w){
  double result=0;
  np::multi_iter iter = np::make_multi_iter(v,w);
  while (iter.not_done())
  {
    double* v_i = reinterpret_cast<double*>(iter.get_data(0));
    double* w_i = reinterpret_cast<double*>(iter.get_data(1));
    result+=(*v_i)*(*w_i);
    iter.next();
  }
  return result;
}


double np_sum(np::ndarray v){
  double result=0;
  np::multi_iter iter = np::make_multi_iter(v);
  while (iter.not_done())
  {
    double* v_i = reinterpret_cast<double*>(iter.get_data(0));
    result+=(*v_i);
    iter.next();
  }
  return result;
}
double np_sum(boost::python::api::object v){
  return np_sum(np::from_object(v));
}

double np_max(boost::python::api::object v){
  double max=0;
  np::multi_iter iter = np::make_multi_iter(np::from_object(v));
  while (iter.not_done())
  {
    double* v_i = reinterpret_cast<double*>(iter.get_data(0));
    max=max<(*v_i)?(*v_i):max;
    iter.next();
  }
  return max;
}

/******/




//#<precipitation related functionality v>


//#NOTE:EPA evaporation in Technical guidance for hazards analysis, eq (7) page G-3
//#TODO:check QS for spilled water
double QR(double u,double BS_s,double T_t){//#in l/day
    double U=u * 1000.0/(3600.0);//#step(wind_speed,t) * 1000.0/(60.0*60.0)#in m/s #km/h->m/s
    double MW=18.01528;//#molecular wight of water in g/mol (gmol often called mol src:https://chemistry.stackexchange.com/questions/53455/g-gmol-vs-lb-lbmol)
    double A=BS_s * 0.00107;//#in ft^2 #cm^2->ft^2
    double VP=pow(10.0,(8.07131-(1730.63/(233.426+T_t-273.15))));//#Vapor pressure by Antoine Equation in mmHg
    double R=82.05; //#in atm cm 3 /g mole
    return ( (0.106 * pow(U,0.78) * pow(MW,(2.0/3.0))* A * VP)/(R* T_t) ) * 453.59237/(1.0/(60.0*24.0)) *1./1000.; //#Rate of release to air ml/day #(Ibs/min) ->  ml/day ->l/day
  }

double QG(double BS_s,double p_t,double t){//#Quantity gathered#in litres
    return (BS_s * p_t*0.1) * 1.0 * 1./1000.;//#*1cm^3=1ml -> l
  }


//      { QG(BS_s,t)-QR(BS_s,T(t))    if 0 < W < BS_c
//  dW= { QG(BS_s,t)                  if W <= 0.0
//      { -QR(BS_s,T(t))               if W >= BS_c
//Note: in the implementation we needed to add functions to make function continuous, otherwise odeint breaks


double dW(double W,double BS_c,double BS_s,double T_t,double p_t,double wss_t,double t){//#in l/day
    double epsilon=1e-3;
    if(0+epsilon < W && W< BS_c-epsilon){
        return QG(BS_s,p_t,t)-QR(wss_t,BS_s,T_t);
      }
    else if(W <= 0.0+epsilon){
        return QG(BS_s,p_t,t) - QR(wss_t,BS_s,T_t)*(W/epsilon);
      }
    else if( W >= BS_c-epsilon){
        return QG(BS_s,p_t,t)*((BS_c-W)/epsilon) - QR(wss_t,BS_s,T_t);
      }
      return 0;//to avoid a warning
}

double a0(double W){
    return 70.0* W;
  }


//#</precipitation related functionality v>


np::ndarray vR_D(double T_t){//#day^-1
    double R=1.987;// # universal gas constant
    np::ndarray vR_D_298K=np::array(p::make_tuple(0.24,0.2088,0.384,0.216,0.372));
    //ro_25_=[0.01066,0.00873,0.01610,0.00898] #replaced by R_D_298K. which:  R_D_298K ~ ro_25*24  #(24 hours)
    np::ndarray vDeltaH_A=np::array(p::make_tuple(10798.0,26018.0,14931.0,15725.0,15725.0));
    np::ndarray vDeltaH_H=np::array(p::make_tuple(100000.0,55990.0,-472379.00,1756481.0,1756481.0));// #-472379 vs. -473379
    np::ndarray vT_1_2H=np::array(p::make_tuple(14184.0,304.6,148.0,447.2,447.2));

    auto ret= vR_D_298K * (T_t/298.0) * np_exp( (vDeltaH_A/R)* ((1.0/298.0)- (1.0/T_t)) ) / ( 1.0+ np_exp( (vDeltaH_H/R)* ((1.0/vT_1_2H)-(1.0/T_t)) ) );
    return np::from_object(ret);
  }


double gamma(double L,double BS,double W){
    double epsilon=1e-4;
    if(BS==0 || W <epsilon)//#W *1000./BS_s <0.1
        return 1.0;//#no water total inhibition#1960 Aedes aegypti (L.), The Yellow Fever Mosquito(Page 165)
    else if(L/BS<=a0(W)-epsilon)
        return 0;
    else if(a0(W)-epsilon < L/BS && L/BS <a0(W)+epsilon){
        //#a (a0-e) + b=0 => b=-a (a0 -e)
        //#a (a0 + e) + b=0.63 => a(a0+e) - a(a0-e) = 2 a e = 0.63 =>a=0.63/(2 e)
        double a=0.63/(2.0*epsilon);
        double b=-a*(a0(W)-epsilon);
        return a * (L/BS) + b;
      }
    else if(L/BS>=a0(W)+epsilon)
        return 0.63;

    return 0;//to avoid a warning
}
np::ndarray vGamma(np::ndarray vL,boost::python::api::object vBS_a,np::ndarray vW){
  //return np.array([gamma(vL[i],vBS_a[i],vW[i]) for i in range(0,len(vL))])
    p::list result=p::list();
    np::multi_iter iter = np::make_multi_iter(vL,np::from_object(vBS_a),vW);
    while (iter.not_done())
    {
      double* l = reinterpret_cast<double*>(iter.get_data(0));
      double* bs_a = reinterpret_cast<double*>(iter.get_data(1));
      double* w = reinterpret_cast<double*>(iter.get_data(2));
      result.append(gamma(*l,*bs_a,*w));
      iter.next();
    }
    return np::array(result);
  }




np::ndarray f(np::ndarray vW, np::ndarray vBS_d){//#TODO:change name to something with meaning
    double epsilon=1e-4;
    np::ndarray vf=np::from_object(vW/(vW+epsilon) * vBS_d);
    if(np_max(vf)<1e-20)
        return vf;
    else
        return np::from_object(vf/np_sum(vf));//#TODO: check this
}

np::ndarray dvE(np::ndarray vE,np::ndarray vL,double A1,double A2,np::ndarray vW,double T_t,double BS_a,np::ndarray vBS_d,double elr,double ovr1,double ovr2){
    double egn=63.0;
    double me=0.01;//#mortality of the egg, for T in [278,303]
    auto ret =egn*( ovr1 *A1  + ovr2* A2)*f(vW,vBS_d) - me * vE - elr* (1-vGamma(vL,BS_a*vBS_d,vW)) * vE;
    return np::from_object(ret);
}

np::ndarray dvL(np::ndarray vE,np::ndarray vL,np::ndarray vW,double T_t,double BS_a,np::ndarray vBS_d,double elr,double lpr){
    double ml=0.01 + 0.9725 * exp(-(T_t-278.0)/2.7035);//#mortality of the larvae, for T in [278,303]
    double alpha0=1.5;//#Parameter to be fitted #1.0#HARDCODED!!!
    double alpha=alpha0/BS_a;//#Σ vα[i] * vL[i]^2= Σ α0/(BS_a* vBS_d[i]) * (L*vBS_d[i])^2 = Σ α0/BS_a * L^2 * vBS_d[i] = α *L^2 *Σ BS_d[i]=α L^2 #Note: on why this is ok even though BS changed.
    auto ret= elr* (1-vGamma(vL,BS_a*vBS_d,vW)) * vE - ml*vL - alpha* vL*vL - lpr *vL ;//#-35.6464*(1.-beta(vW,vBS_od,vBS_id))*L#-24.*(1.-beta(vW))*L# -log(1e-4/5502.)/(1.)=17.823207313460703
    return np::from_object(ret);
}

np::ndarray dvP(np::ndarray vL, np::ndarray vP, double T_t, double lpr, double par){
    double mp=0.01 + 0.9725 * exp(-(T_t-278.0)/2.7035);//#death of pupae
    auto ret= lpr*vL - mp*vP  - par*vP;
    return np::from_object(ret);
}

double dA1(np::ndarray vP, double A1, double T_t, double par, double ovr1){
    double ef=0.83;//#emergence factor
    double ma=0.09;//#for T in [278,303]
    return np_sum(np::from_object(par*ef*vP/2.0)) - ma*A1 - ovr1*A1;
}

double dA2(double A1, double A2, double T_t, double ovr1){
    double ma=0.09;//#for T in [278,303]
    return ovr1*A1 - ma*A2;
}

void initialize(){
  Py_Initialize();
  np::initialize();
}

/************************************************/
//this main method has no effect on python, is just to test from a c++ environment


#endif
int main(){

  std::cout << dA2(1,2,3,4);
}
