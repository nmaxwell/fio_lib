#ifndef BFIO_INTERPOLATION_H
#define BFIO_INTERPOLATION_H

#include <math.h>
#include <stdlib.h>

double lagrange_basis(double point, int size_support, double *support, int index);

double lagrange_interpolate(double point, int size_support, double *support, double *function_values);

int lagrange_matrix(double *matrix, int n_input_points, int n_output_points, double *input_points, double *output_points);

int chebyshev_grid(double *grid, int n_points, double start, double end);


#endif
