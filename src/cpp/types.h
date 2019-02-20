#ifndef TypesH
#define TypesH

#include <valarray>

typedef double scalar;
typedef std::valarray<scalar> tensor;//TODO:find a better name.This is a vector, with usual scalar-vector product and vector-vector sum, with addition to a vector-vector pointwise product.
typedef std::valarray<tensor> matrix;
matrix operator*(scalar lhs,const matrix& rhs){
    unsigned int m=rhs.size();
    unsigned int n=rhs[0].size();
    matrix M=matrix(tensor(n),m);
    for(unsigned int i=0;i<m;i++) M[i]=lhs*rhs[i];
    return M; // return the result by value (uses move constructor)
  }
#endif
