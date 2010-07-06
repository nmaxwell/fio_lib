#ifndef EXACT_FIO_CPP
#define EXACT_FIO_CPP

#include "exact_fio.h"


#define FIO_CONST_2PI 6.2831853071795865

int exact_fio_2d( double complex *output, int N, double complex *input, double (*phase)(double, double, double, double) )
{
    if (output==NULL) return 1;
    if (N==0) return 2;
    
	double halfN = ((double)N)/2;
    
    for (int i1=0; i1<N; i1++)
	for (int i2=0; i2<N; i2++)
	{
		double complex sum=0.0*I;
			
		for (int j1=0; j1<N; j1++)
		for (int j2=0; j2<N; j2++)
			sum += input[j2+j1*N]*cexp(I*FIO_CONST_2PI*( phase((double)i1/N, (double)i2/N, (double)j1-halfN, (double)j2-halfN) ));
	
		output[i2+i1*N] = sum;
	}
    
    return 0;
}


int exact_fio_dft( double complex *output, int N, double complex *input )
{
    if (output==NULL) return 1;
    if (N==0) return 2;
    
    for (int i1=0; i1<N; i1++)
    for (int i2=0; i2<N; i2++)
    {
        double complex sum=0.0*I;
        
        for (int j1=0; j1<N; j1++)
        for (int j2=0; j2<N; j2++)
            sum += input[j2+j1*N]*cexp(I*FIO_CONST_2PI*( -((double)(i1*j1+i2*j2))/N ));
        
        output[i2+i1*N] = sum;
    }
    
    return 0;
}



#endif
