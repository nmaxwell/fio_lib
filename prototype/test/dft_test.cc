
#include <time.h>
#include <math.h>
#include <sys/stat.h>

#include "../c_code/exact_fio.h"

//#include "fftw3.h"

#define _2pi 6.2831853071795865


double get_real_time()
{
	return (double)( clock() )/(CLOCKS_PER_SEC);
}



int N=32;


double dft_phase(double x1, double x2, double k1, double k2)
{
	double halfN = ((double)N)/2;
	
    return -(x1*(k1+halfN) + x2*(k2+halfN));
}


double uniform(double min, double max)
{
    return (max-min)*((double)rand())/RAND_MAX+min;
}

int main()
{
	//printf("%d\n",  CLOCKS_PER_SEC );
    srand(time(NULL));
    
<<<<<<< HEAD
    N=64;
=======
    N=32;
>>>>>>> 3986310ff80c68cafd3241fa9b449d03dfcc2160
    
    complex double *input= (complex double *)malloc(N*N*sizeof(complex double));
    complex double *temp1= (complex double *)malloc(N*N*sizeof(complex double));
    complex double *temp2= (complex double *)malloc(N*N*sizeof(complex double));
    //complex double *fftw_output= (complex double *)malloc(N*N*sizeof(complex double));
	
    complex double *fio_exact_output= (complex double *)malloc(N*N*sizeof(complex double));
    complex double *dft_output= (complex double *)malloc(N*N*sizeof(complex double));
    
    for (int k=0; k<N*N; k++)
    {
        //fftw_output[k]=0;
        fio_exact_output[k]=0;
        
        input[k] = I*uniform(-1.0,1.0)+uniform(-1.0,1.0);
    }
    
    
    int err=0;
	double t1,t2, time_exact_fio_2d, time_exact_fio_dft;
	
	t1 = get_real_time();
    err += exact_fio_2d( fio_exact_output, N, input, dft_phase );
    t2 = get_real_time();
	time_exact_fio_2d = t2-t1;
	
	
    //fftw_plan plan = fftw_plan_dft_2d( N, N, temp1, temp2, FFTW_FORWARD, FFTW_ESTIMATE );
    //fftw_execute_dft( plan, input, fftw_output );
    
	t1 = get_real_time();
    err += exact_fio_dft( dft_output, N, input );
    t2 = get_real_time();
	time_exact_fio_dft = t2-t1;
	
	printf("time for exact_fio_2d:\t%f\ntime for exact_fio_dft:\t%f\n", time_exact_fio_2d, time_exact_fio_dft);
	
	if (err != 0) printf("error\n");
    
    double EX2=0;
    double EX=0;
    double min,max;
    
    double magnitude=0;
    for (int k=0; k<N*N; k++)
        magnitude += cabs(dft_output[k])*cabs(dft_output[k])/(N*N);
    magnitude = sqrt(magnitude);
	
    for (int k=0; k<N*N; k++)
    {
        //double x = log10( fabs(creal(fftw_output[k])-creal(fio_output[k])) + fabs(cimag(fftw_output[k])-cimag(fio_output[k]))  +1E-17 );
        
        double x = log10( fabs(cabs(dft_output[k] - fio_exact_output[k]))/magnitude +1E-17 );
        
        if (k==0)
        {
            min=x;
            max=x;
        }
        else
        {
            if (min>x) min=x;
            if (max<x) max=x;
        }
        
        //printf( "%f\n", x  );
        
        EX += x/(N*N);
        EX2 += x*x/(N*N);
    }
    
    printf( "mean: %f\t var: %f\t min: %f\t max: %f\n", EX, EX2-EX*EX, min, max  );
    
    
    
    //fftw_destroy_plan( plan );
    free(input);
    free(temp1);
    free(temp2);
    //free(fftw_output);
    free(fio_exact_output);
    
    return 0;
}




