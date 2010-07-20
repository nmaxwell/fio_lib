
//#include <google/profiler.h>

#if defined(__APPLE__)
#include <CoreServices/CoreServices.h>
#endif


#include "../cppbfio/bfio_prototype_1.cpp"
//#include "../cppbfio/bfio_lexing.cpp"
#include "../cppbfio/exact_fio.cpp"
#include <iostream>
#include <time.h>
#include <stdio.h>
#include <fstream>

using namespace std;


double get_real_time()
{

#if defined(__APPLE__)
 
  Nanoseconds Nsecs = AbsoluteToNanoseconds(UpTime());
  return (double)(UnsignedWideToUInt64(Nsecs))*1.0E-9;

#elif defined(__linux)

  timespec current_time;
  clock_gettime(CLOCK_REALTIME, &current_time);
  return (double) ((double)current_time.tv_sec+(double)current_time.tv_nsec/(1.0E9));
#else

  return 0;

#endif

}

double uniform(double low, double high)
{
  return (high-low)*(double)(rand()%(RAND_MAX))/(RAND_MAX-1)+low;
}



int N;


double phase(double x1, double x2, double k1, double k2)
{
  return (x1*k1+x2*k2);



  double norm_k = sqrt(k1*k1+k2*k2);
  double angle = atan2(k2,k1);
  // f(x1,x2,angle)
  
  double f = fabs(x1)*x1*sin(x2*6)*exp(cos(angle*3));
  
  return norm_k*f;
}


double phase_odd(double x1, double x2, double k1, double k2)
{
  double halfN = 0.5*N;
  
  return phase(x1/N, x2/N, k1-halfN, k2-halfN);
}

complex<double> amplitude(double x1, double x2, double k1, double k2)
{
  return 1.0;
  return exp(-(x1*x2+k1*k2)*10.0);
}

complex<double> amplitude_odd(double x1, double x2, double k1, double k2)
{
  double halfN = 0.5*N;
  
  return amplitude(x1/N, x2/N, k1-halfN, k2-halfN);
}



int test0(const char *input_fname )
{
  int error=0;
  
  srand ( time(NULL) );

  ifstream input_file( input_fname );
  int m,n;
  input_file.read((char*)(&m), 4);
  input_file.read((char*)(&n), 4);
  assert( n == m );
  N = m;
  cout << "N: " << N << endl;
  complex<double>  *input = (complex<double> *)calloc(N*N, 16);
  complex<double> *output = (complex<double> *)calloc(N*N, 16);

  input_file.read((char *)(input), m*n*16);
  
  /*
    int x1 = 1;
    int x2 = 2;
  
    for (int k1=0; k1<N; k1++)
    for (int k2=0; k2<N; k2++)
    input[k1+k2*N] = exp(FIO_CONST_2PI_IMAG*(double)(x1*k1+x2*k2)/(double)N);
  */
  
  cout << endl;
  //error = exact_dfio_2d(  N, output, input, phase );
  
  for (int i1=0; i1<4; i1++)
    for (int i2=0; i2<4; i2++)
      cout << i1 << "\t" << i2 << "\t" << output[i1+i2*N]/(double)(N*N) << endl;
  
  
  
  free(input);
  free(output);
  return error;
}



int test1(const char *input_fname )
{
  int error=0;
  
  srand ( time(NULL) );

  ifstream input_file( input_fname );
  int m,n;
  input_file.read((char*)(&m), 4);
  input_file.read((char*)(&n), 4);
  assert( n == m );
  N = m;
  cout << "N: " << N << endl;
  complex<double>  *input = (complex<double> *)calloc(N*N, 16);
  complex<double> *output = (complex<double> *)calloc(N*N, 16);

  input_file.read((char *)(input), m*n*16);  
  
  cout << endl << endl;
  int log2N = (int)(log(N)/log(2));
  int start_level = log2N-3;
  int end_level = 3;
  int n_cheby = 5;

  //  ProfilerStart("/Users/maxwelna/fio_lib/prototype/test/profile_out");
  double start_time = get_real_time();

  error = bfio_core(input, output, N, start_level, end_level, n_cheby, phase_odd, amplitude_odd);
  if (error) return error;
  
  double end_time = get_real_time();
  // ProfilerStop();
  double bfio_run_time = end_time - start_time;
  
  
  
  int n_samples = 1 << 5;
  int *sample_indexes = new int [n_samples*2];
  complex<double> *exact_output = (complex<double> *)calloc(n_samples, 16);
  for (int k=0; k<n_samples; k++) {
    sample_indexes[2*k] = rand()%N;
    sample_indexes[2*k+1] = rand()%N; }
  
  start_time = get_real_time();
  error = exact_dfio_2d_sample(n_samples, sample_indexes, N,  exact_output, input, phase, amplitude);
  if (error) return error;
  end_time = get_real_time();

  double fio_run_time = (end_time-start_time)*(N*N)/(double)(n_samples);
  double err=0.0;
  for (int k=0; k<n_samples; k++) {
    err += pow(abs((output[sample_indexes[2*k]+N*sample_indexes[2*k+1]]-exact_output[k])),2);
        cout << k << "\t" << sample_indexes[k*2] << "\t" << sample_indexes[2*k+1] << "\t" << output[sample_indexes[2*k]+N*sample_indexes[2*k+1]] << "\t" << exact_output[k] << endl;
 }
  double mag=0.0;
  for (int k=0; k<n_samples; k++)
    mag += pow(abs(exact_output[k]),2);
  err = sqrt(err/mag);

  cout << "BFIO run time: " << bfio_run_time  << endl;
  cout << "FIO run time: " << fio_run_time << endl;
  cout << "Speedup: " << fio_run_time/bfio_run_time << endl;
  cout << "Relative l2 error: " << err << endl;
  
  
  delete [] sample_indexes;
  free(input);
  free(output);
  return error;
}






int main()
{
  int error = 0;

    //error = test0();


    
  //while (1)
  error = test1("../lexing/input/f_256.bin");

  if (error)
    printf("error ocurred: %d\n", error);
  else
    printf("All good.\n");
    
    
  return 0;
}
