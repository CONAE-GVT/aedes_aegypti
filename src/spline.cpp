//>g++ -Wall -std=c++11 src/spline.cpp
#include <iostream>
#include <vector>
#include "spline.h"

int main(int argc, char** argv) {

   std::vector<double> X(500), Y(500);
   for(int i=0;i<500;i++){
      X[i]=float(i);
      Y[i]=X[i]*X[i];
   }
   tk::spline s=tk::spline();
   s.set_points(X,Y);    // currently it is required that X is already sorted

   double x=1.5;
   std::cout << "spline at "<< x<< " is "<<s(x)<<std::endl;

   return EXIT_SUCCESS;
}
