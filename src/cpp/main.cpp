//g++ -Wall main.cpp
//>./a.out > lorenz.dat
//gnuplot>plot "lorenz.dat" using 1:2 with lines
//https://www.geeksforgeeks.org/std-valarray-class-c/
#include <valarray>
#include <iostream>
#include "weather.h"
#include "rk.h"
#include "equations.h"

int main(){
  /*
  //vectors
  std::valarray<double> vA = { 10, 2, 20, 1, 30 };
  std::valarray<double> vB = { 11, 2, 20, 1, 30 };
  auto vC=vA*vB + vA + (2.*vB);
  std::cout <<"vC[0]= "<<vC[0]<<std::endl;

  //matrix
  std::valarray<std::valarray<double>> mA = {{ 10, 2, 20, 1, 30 },{ 10, 2, 20, 1, 30 }};
  std::valarray<std::valarray<double>> mB=mA*mA;
  std::cout << "mB[0,0]=(mA*mA)[0,0]= "<<mB[0][0]<<std::endl;
  //matrix muting
  mA[0]={ 11, 2, 20, 1, 30 };
  std::cout<<"mA= "<<std::endl;
  for(auto row:mA){
    for(auto e:row){
      std::cout << e<<" ";
    }
    std::cout << std::endl;
  }
  mB=mA*mA;
  std::cout << "mB[0,0]= "<<mB[0][0]<<std::endl;
  */

  Weather weather=Weather("data/public/wunderground.csv", "2015-7-1", "2018-7-1");
  //std::cout << weather.T(2.5)<<std::endl;

  //TODO:HARDCODED! find a better way of doing a timerange
  std::vector<double> time_range;
  double h=1/40.;
  for(unsigned int i=0;i<60/h;i++) time_range.push_back(i*h);

  std::valarray<double> Y0={ 10.0 , 1.0 , 1.0 };
  std::vector<std::valarray<double>> Y=RK::solve(diff_eqs,Y0,time_range,20);
  for(unsigned int i=0;i<time_range.size();i++) std::cout << time_range[i] << '\t' << Y[i][0] << '\t' << Y[i][1] << '\t' << Y[i][2] << std::endl;
}
