//g++ -Wall main.cpp
//>./a.out > lorenz.dat
//gnuplot>plot "lorenz.dat" using 1:2 with lines
//https://www.geeksforgeeks.org/std-valarray-class-c/
#include<vector>
#include <valarray>
#include <iostream>
#include "otero_precipitation.h"
int main(){
  Model model=Model();
  std::vector<std::valarray<double>> Y=model.solveEquations();
  std::vector<double> time_range=model.time_range;
  for(unsigned int i=0;i<time_range.size();i++) std::cout << time_range[i] << '\t' << Y[i][0] << '\t' << Y[i][1] << '\t' << Y[i][2] << std::endl;
}
