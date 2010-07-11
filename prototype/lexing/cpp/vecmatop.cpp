#include "blas.h"
#include "lapack.h"
#include "numvec.hpp"
#include "nummat.hpp"
#include "numtns.hpp"

#include "vecmatop.hpp"

using std::cerr;


//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C)
{
  assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
  char transa = 'N';
  char transb = 'N';
  int m = C.m();  int n = C.n();  int k = A.n();
  dgemm_(&transa, &transb, &m, &n, &k,
	 &alpha, A.data(), &m, B.data(), &k, &beta, C.data(), &m);
  return 0;
}

// ---------------------------------------------------------------------- 
int zgemm(cpx alpha, const CpxNumMat& A, const CpxNumMat& B, cpx beta, CpxNumMat& C)
{
  assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
  char transa = 'N';
  char transb = 'N';
  int m = C.m();  int n = C.n();  int k = A.n();
  zgemm_(&transa, &transb, &m, &n, &k,
	 &alpha, A.data(), &m, B.data(), &k, &beta, C.data(), &m);
  return 0;
}

//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y)
{
  assert(Y.m() == A.m());
  assert(A.n() == X.m());
  char trans = 'N';
  int m = A.m();  int n = A.n();
  int incx = 1;  int incy = 1;
  dgemv_(&trans, &m, &n, &alpha, A.data(), &m, X.data(), &incx, &beta, Y.data(), &incy);
  return 0;
}

//Y <- a M X + b Y
// ----------------------------------------------------------------------
int zgemv(cpx alpha, const CpxNumMat& A, const CpxNumVec& X, cpx beta, CpxNumVec& Y)
{
  assert(Y.m() == A.m());
  assert(A.n() == X.m());
  char trans = 'N';
  int m = A.m();  int n = A.n();
  int incx = 1;  int incy = 1;
  zgemv_(&trans, &m, &n, &alpha, A.data(), &m, X.data(), &incx, &beta, Y.data(), &incy);
  return 0;
}

//-------------------------------------------------------------------------
//just the transpose,
int ztran(const CpxNumMat& A, CpxNumMat& B)
{
  B.resize(A.n(), A.m());
  for(int i=0; i<B.m(); i++)
    for(int j=0; j<B.n(); j++)
      B(i,j) = A(j,i);
  return 0;
}
