//g++ -Wall main.cpp
//>./a.out > lorenz.dat
//gnuplot>plot "lorenz.dat" using 1:2 with lines
//https://www.geeksforgeeks.org/std-valarray-class-c/
//cppcheck --enable=all src/cpp/main.cpp
#include<vector>
#include <iostream>
#include "types.h"
#include "otero_precipitation.h"

int main(){
  Model model=Model();
  std::vector<tensor> Y=model.solveEquations();
  std::vector<scalar> time_range=model.time_range;
  for(unsigned int i=0;i<time_range.size();i++){
    std::cout << time_range[i] << '\t';
    for(unsigned int j=0;j<Y[i].size();j++) std::cout << Y[i][j] << '\t';
    std::cout <<  std::endl;
  }
  //W
  /*
  for(unsigned int i=0;i<time_range.size();i++){
    std::cout << time_range[i] << '\t';
    tensor vW_t=model.parameters.vW(time_range[i]);
    for(unsigned int j=0;j<vW_t.size();j++) std::cout << vW_t[j] << '\t';
    std::cout <<  std::endl;
  }
  */
}
