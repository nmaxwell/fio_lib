#ifndef EXACT_FIO_H
#define EXACT_FIO_H

/*

documentation:

    nicholas.maxwell@gmail.com for questions.
    
    This file is for evaluating the Fourier Integral Operator
    
    u(x) = \sum_{k \in \Omega} f(k)  e^{ 2 \pi i \Phi(x,k)}, x \in X
    
    where i*i=-1, and \Phi is homogenous in the second variable, i.e.
        
        \Phi(x,\lambda k) = \lambda \Phi(x,k), for all \lambda > 0
    
    and X = \{ (i_1/N, i_2/N) ;  0 \le i_1, i_2 < N \},
    \Omega = \{ (k_1, k_2) ;  -N/2 \le k_1, k_2 < N/2 \}
    
    for some positive integer N, here a power of two.
    
    naming:
        
        the function 'f' in the above formula is refered to as input,
        and the function \Phi is refered to as the phase function.
    
    arrays are accessed as row-major, c-style
*/

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


<<<<<<< HEAD:prototype/cbfio/exact_fio.h
int exact_fio_2d(double complex *output, int N, double complex *input, double (*phase)(double, double, double, double) );

int exact_fio_dft(double complex *output, int N, double complex *input);
=======
//using namespace std;

int exact_fio_2d(double complex *output, int N, double complex *input, double (*phase)(double, double, double, double) );


int exact_fio_dft(double complex *output, int N, double complex *input);


>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160:prototype/c_code/exact_fio.h





#endif
