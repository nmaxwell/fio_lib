#ifndef BFIO_INTERPOLATION_H
#define BFIO_INTERPOLATION_H

#include <math.h>


double lagrange_basis(double point, int size_support, double *support, int index);

double lagrange_interpolate(double point, int size_support, double *support, double *function_values);

double * chebyshev_grid(int n_points, double start, double end);




#endif
