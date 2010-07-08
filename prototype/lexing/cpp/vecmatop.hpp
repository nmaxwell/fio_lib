#ifndef _VECMATOP_HPP_
#define _VECMATOP_HPP_

#include "nummat.hpp"

//--------------------------------------------------
int dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C);

int dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y);

//--------------------------------------------------
int zgemm(cpx alpha, const CpxNumMat& A, const CpxNumMat& B, cpx beta, CpxNumMat& C);

int zgemv(cpx alpha, const CpxNumMat& A, const CpxNumVec& X, cpx beta, CpxNumVec& Y);

int ztran(const CpxNumMat& A, CpxNumMat& B);


#endif
