#ifndef BFIO_PROTOTYPE_H
#define BFIO_PROTOTYPE_H

#define _debug_here( pos ) ( printf( "debug: file %s; \n line %d; \t code %d\n", __FILE__ , __LINE__, pos ) );

/*
 The goal here is transparecy and reasonable efficiency.
 
 
 
 TODO:
    Set up error codes and handling, etc.
 
 
*/

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bfio_interpolation.h"
    
#define bfio_pi 3.1415926535897932384626433832795028841
#define bfio_2pi 6.2831853071795865
#define bfio_pisqrt2 4.4428829381583661
#define bfio_sqrt2_2 0.70710678118654757


int bfio_pow2(int k)
{
    return 1<<k;
}


// typedefs

typedef struct {double x1; double x2; double complex value;} bfio_data_point;

typedef struct {int size; bfio_data_point *data;} bfio_data_list;

typedef struct
{
    int N;
    int n_cheby;
    int start_level;
    int end_level;
    double (*phase)(double, double, double, double);
    

    bfio_data_list *input_data;
    int input_data_size;

    
    double *unit_cheby_grid;
    
    
} bfio_session;




// functions in earnest

int bfio_initialize_session(bfio_session *session,  int N, int n_cheby, int start_level, int end_level, double (*phase)(double, double, double, double));

int bfio_destruct_session(bfio_session *session);


int bfio_initialize_input_data(bfio_session *session, int j1, int j2, int size);

int bfio_free_input_data(bfio_session *session, int j1, int j2);

int bfio_set_input_data_array(bfio_session *session, int j1, int j2, int size, double *x1, double *x2, double complex *value);

int bfio_set_input_data_point(bfio_session *session, int j1, int j2, int position, double x1, double x2, double complex value);


int bfio(bfio_session *session);






#endif
