#ifndef EXACT_FIO_CPP
#define EXACT_FIO_CPP

#include "exact_fio.h"


#define FIO_CONST_2PI 6.2831853071795865
#define FIO_CONST_2PI_IMAG (complex<double> (0., FIO_CONST_2PI))

#ifdef __cplusplus
extern "C" {
#endif

  int exact_dfio_2d_ca(  int N, complex<double> *output, complex<double> *input, double (*phase)(double, double, double, double) )
  {
    if (output==NULL) return 1;
    if (input==NULL) return 1;
    if (N==0) return 2;
    
    double halfN = ((double)N)/2;
    
    for (int i1=0; i1<N; i1++)
      for (int i2=0; i2<N; i2++)    {
        complex<double> sum=0.0;
            
        for (int j1=0; j1<N; j1++)
	  for (int j2=0; j2<N; j2++)
            sum += input[j1+j2*N]*exp(FIO_CONST_2PI_IMAG*( phase((double)i1/N, (double)i2/N, (double)j1-halfN, (double)j2-halfN) ));
    
        output[i1+i2*N] = sum;
      }
    
    return 0;
  }




  int exact_dfio_2d_ca_sample( int n_samples, int *sample_indexes, int N, complex<double> *output, complex<double> *input, double (*phase)(double, double, double, double) )
  {
    if (output==NULL) return 1;
    if (input==NULL) return 1;
    if (sample_indexes==NULL) return 1;
    if (N==0) return 2;
    
    double halfN = ((double)N)/2;
    

    for (int k=0; k<n_samples; k++) {
      int i1,i2;
      i1 = sample_indexes[2*k];
      i2 = sample_indexes[2*k+1];
      
      complex<double> sum=0.0;
      
      for (int j1=0; j1<N; j1++)
	for (int j2=0; j2<N; j2++)
	  sum += input[j2+j1*N]*exp(FIO_CONST_2PI_IMAG*( phase((double)i1/N, (double)i2/N, (double)j1-halfN, (double)j2-halfN) ));
    
      output[k] = sum;
    }
    
    return 0;
  }



  int exact_dft_2d( int N, complex<double> *output, complex<double> *input )
  {
    if (output==NULL) return 1;
    if (input==NULL) return 1;
    if (N==0) return 2;
    
    for (int i1=0; i1<N; i1++)
      for (int i2=0; i2<N; i2++)    {
        complex<double> sum=0.0;
        
        for (int j1=0; j1<N; j1++)
	  for (int j2=0; j2<N; j2++)
            sum += input[j1+j2*N]*exp(FIO_CONST_2PI_IMAG*( -((double)(i1*j1+i2*j2))/N ));
        
        output[i1+i2*N] = sum;
      }
    
    return 0;
  }


  int exact_dft_2d_sample( int n_samples, int *sample_indexes, int N, complex<double> *output, complex<double> *input)
  {
    if (output==NULL) return 1;
    if (input==NULL) return 1;
    if (sample_indexes==NULL) return 1;
    if (N==0) return 2;
    
    for (int k=0; k<n_samples; k++) {
      int i1,i2;
      i1 = sample_indexes[2*k];
      i2 = sample_indexes[2*k+1];
      
      complex<double> sum=0.0;
      
      for (int j1=0; j1<N; j1++)
	for (int j2=0; j2<N; j2++)
	  sum += input[j1+j2*N]*exp(FIO_CONST_2PI_IMAG*( -((double)(i1*j1+i2*j2))/N ));
    
      output[k] = sum;
    }
    
    return 0;
  }





  int exact_dfio_2d(  int N, complex<double> *output, complex<double> *input, double (*phase)(double, double, double, double), complex<double> (*amplitude)(double, double, double, double))
  {
    if (output==NULL) return 1;
    if (input==NULL) return 1;
    if (N==0) return 2;
    
    double halfN = ((double)N)/2;
    
    for (int i1=0; i1<N; i1++)
      for (int i2=0; i2<N; i2++)    {
        complex<double> sum=0.0;
            
        for (int j1=0; j1<N; j1++)
	  for (int j2=0; j2<N; j2++)
            sum += input[j1+j2*N]*exp(FIO_CONST_2PI_IMAG*( phase((double)i1/N, (double)i2/N, (double)j1-halfN, (double)j2-halfN) ))*amplitude( (double)i1/N, (double)i2/N, (double)j1-halfN, (double)j2-halfN );
    
        output[i1+i2*N] = sum;
      }
    
    return 0;
  }




  int exact_dfio_2d_sample( int n_samples, int *sample_indexes, int N, complex<double> *output, complex<double> *input, double (*phase)(double, double, double, double), complex<double> (*amplitude)(double, double, double, double) )
  {
    if (output==NULL) return 1;
    if (input==NULL) return 1;
    if (sample_indexes==NULL) return 1;
    if (N==0) return 2;
    
    double halfN = ((double)N)/2;

    for (int k=0; k<n_samples; k++) {
      int i1,i2;
      i1 = sample_indexes[2*k];
      i2 = sample_indexes[2*k+1];
      
      complex<double> sum=0.0;
      
      for (int j1=0; j1<N; j1++)
	for (int j2=0; j2<N; j2++)
	  sum += input[j2+j1*N]*exp(FIO_CONST_2PI_IMAG*( phase((double)i1/N, (double)i2/N, (double)j1-halfN, (double)j2-halfN) ))*amplitude( (double)i1/N, (double)i2/N, (double)j1-halfN, (double)j2-halfN );
    
      output[k] = sum;
    }
    
    return 0;
  }















#ifdef __cplusplus
}
#endif

#endif
