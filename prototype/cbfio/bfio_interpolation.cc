#ifndef BFIO_INTERPOLATION_CC
#define BFIO_INTERPOLATION_CC

#include "bfio_interpolation.h"

double lagrange_basis(double point, int size_support, double *support, int index)
{
    if (support==0) return 0;
    
    double prod=1.0;
    for (int k=0; k<size_support; k++)
        if (k != index)
            prod *= (point - support[k])/(support[index]-support[k]);
    
    return prod;
}

double lagrange_interpolate(double point, int size_support, double *support, double *function_values)
{
    if (support==0 || function_values==0) return 0;
    
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

double * chebyshev_grid(int n_points, double start, double end)
{
    double *ptr = (double *)calloc(n_points, sizeof(double));
    if (ptr == NULL)
        return ptr;
    
    for (int k=0; k<n_points; k++)
        ptr[n_points-k-1] = (cos((double)k*3.1415926535897932384/(n_points-1))*0.5+0.5)*(end-start)+start;
    
    return ptr;
}



#endif
