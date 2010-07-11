#ifndef BFIO_INTERPOLATION_CPP
#define BFIO_INTERPOLATION_CPP

#include "bfio_interpolation.h"

double lagrange_basis(double point, int size_support, double *support, int index)
{
    if (support==NULL) return 0;
    
    double prod=1.0;
    for (int k=0; k<size_support; k++)
        if (k != index)
            prod *= (point - support[k])/(support[index]-support[k]);
    
    return prod;
}

double lagrange_interpolate(double point, int size_support, double *support, double *function_values)
{
    if (support==NULL || function_values==NULL) return 0;
    
    double sum=0;
    
    for (int j=0; j<size_support; j++)
    {
        double prod=1.0;
        for (int k=0; k<size_support; k++)
            if (k != j)
                prod *= (point - support[k])/(support[j]-support[k]);
        
        sum += prod*function_values[j];
    }
    return sum;
}

int lagrange_matrix(double *matrix, int n_input_points, int n_output_points, double *input_points, double *output_points)
{
  if (matrix == NULL) return 1;
  
  for (int j=0; j<n_input_points; j++)
    for (int i=0; i<n_output_points; i++)
      matrix[i*n_input_points+j] = lagrange_basis(output_points[i], n_input_points, input_points, j);

  return 0;
}



int chebyshev_grid(double *grid, int n_points, double start, double end)
{
  if (grid == NULL) return 1;

  for (int k=0; k<n_points; k++)
    grid[n_points-k-1] = (cos((double)k*3.1415926535897932384/(n_points-1))*0.5+0.5)*(end-start)+start;
    
  return 0;
}







#endif
