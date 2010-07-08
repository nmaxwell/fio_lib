#ifndef BFIO_PROTOTYPE_H
#define BFIO_PROTOTYPE_H

#define _debug_here( pos ) ( printf( "debug: file %s; \n line %d; \t code %d\n", __FILE__ , __LINE__, pos ) );

/*
 The goal here is transparecy and reasonable efficiency.
 
 
 
 TODO:
    Set up error codes and handling, etc.
 
 
*/

#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
#include <assert.h>

using namespace std;

#include "bfio_interpolation.cpp"
    
#define bfio_pi 3.1415926535897932384626433832795028841
#define bfio_2pi 6.2831853071795865
#define bfio_pisqrt2 4.4428829381583661
#define bfio_sqrt2_2 0.70710678118654757

#define bfio_iu complex<double> (0,1)


int bfio_pow2(int k)
{
    return 1<<k;
}


class bfioDataPoint
{
 public:
  double x1;
  double x2;
  complex<double> value;

 bfioDataPoint(double const &x1, double const &x2, complex<double> const &value):x1(x1), x2(x2), value(value) {};
 bfioDataPoint():x1(0),x2(0),value(0) {};
 bfioDataPoint(const bfioDataPoint &rhs):x1(rhs.x1), x2(rhs.x2), value(rhs.value) {};

};


class bfioSession
{
 public:
  int N;
  int n_cheby;
  int start_level;
  int end_level;
  double (*phase)(double, double, double, double);

  double *unit_cheby_grid;

 private:
  vector<bfioDataPoint> *input_data;
  int input_data_size;

 public:
  bfioSession();
  bfioSession(const bfioSession &);
  bfioSession(int N, int n_cheby, int start_level, int end_level, double (*phase)(double, double, double ,double));
  ~bfioSession();

 private:
  int initialize();
  int bfioSession_free();

 public:
  vector<bfioDataPoint>& input(int index1, int index2);
  int execute();
};


int bfio_core(complex<double> *coefficients, int N, int start_level, int end_level, int n_cheb, double (*phase)(double, double, double ,double));





#endif
