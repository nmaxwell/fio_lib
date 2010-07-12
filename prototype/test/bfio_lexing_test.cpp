
#include "../cppbfio/bfio_lexing.cpp"
#include <iostream>
#include <time.h>
#include <stdio.h>

double get_real_time()
{
  /*timespec current_time;
	clock_gettime(CLOCK_REALTIME, &current_time);
	return (double) ((double)current_time.tv_sec+(double)current_time.tv_nsec/(1.0E9)); */
  return clock();
}


int N;

double dft_phase(double x1, double x2, double k1, double k2)
{
  double halfN = ((double)N)/2;
	
  return -(x1*(k1+halfN) + x2*(k2+halfN));
}


double uniform(double low, double high)
{
  return (high-low)*(double)(rand()%(RAND_MAX))/(RAND_MAX-1)+low;
}



int test1(char input_fname[] )
{
  int error=0;
  
  srand ( time(NULL) );

  ifstream input_file;
  int m,n;
  input_file.read((char*)(&m), 4);
  input_file.read((char*)(&n), 4);
  assert( n == m );
  N = m;

  complex<double>  *input = (complex<double> *)calloc(N*N, 16);
  complex<double> *output = (complex<double> *)calloc(N*N, 16);

  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      input_file.read(input+j*N+i, 16);
  
  int log2N = (int)(log(N)/log(2));
  int start_level = log2N=3;
  int end_level = 3;
  int n_cheby = 9;

  error = bfio_lexing(input, output, N, start_level, end_level, n_cheby, dft_phase);
  

  free(input);
  free(output);
  return error;
}






int main()
{
  int error = 0;
    
  //while (1)
  error = test1();
    
  if (error)
    printf("error ocurred: %d\n", error);
  else
    printf("All good.\n");
    
    
  return 0;
}
