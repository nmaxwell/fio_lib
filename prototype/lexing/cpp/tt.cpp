#include "bfio.hpp"
#include "serialize.hpp"

#include <CoreServices/CoreServices.h>

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





int optionsCreate(int argc, char** argv, map<string,string>& options)
{
  options.clear();
  //2. get extra data
  for(int k=1; k<argc; k=k+2) {
    options[ string(argv[k]) ] = string(argv[k+1]);
  }
  return 0;
}

double phase(double x1, double x2, double k1, double k2)
{
  double norm_k = sqrt(k1*k1+k2*k2);
  double angle = atan2(k2,k1);
  // f(x1,x2,angle)
  
  double f = fabs(x1)*x1*sin(x2*6)*exp(cos(angle*3));
  
  return norm_k*f;
}


int main(int argc, char** argv)
{
  //0. init and get options
  srand48(time(NULL));
  clock_t ck0, ck1;
  time_t t0, t1;
  map<string,string> opts;
  optionsCreate(argc, argv, opts);
  //1. data
  vector<int> all(1,1);
  map<string,string>::iterator mi;
  int N;
  mi = opts.find("-N");
  if(mi!=opts.end()) {    istringstream ss((*mi).second);    ss>>N;  }
  mi = opts.find("-ffile");
  char ffile[100];  {istringstream ss((*mi).second);	ss>>ffile;}
  CpxNumMat f;  {    ifstream fin(ffile);    iC( deserialize(f, fin, all) );  }
  CpxNumMat u(N,N);  setvalue(u,cpx(0,0));
  //2. bfio
  BFIO bfio("bfio_");
  bfio.test_phase = phase;
  iC( bfio.setup(opts) );
  //
  double start_time = get_real_time();
  double time_eval;
  if(N<=256) {
    ck0 = clock();
    iC( bfio.eval(f,u) );
    ck1 = clock();    time_eval = double(ck1-ck0)/CLOCKS_PER_SEC;
  } else {
    t0 = time(0);
    iC( bfio.eval(f,u) );
    t1 = time(0);    time_eval = difftime(t1,t0);
  }
  double end_time = get_real_time();
  /*
  double relerr = 0;
  int NC = 64;
  ck0 = clock();
  iC( bfio.check(f,u,NC,relerr) );
  ck1 = clock();
  double time_chck = double(ck1-ck0)/CLOCKS_PER_SEC * N * N / double(NC);
  */
  //printf("RESULT\n");
  printf("N  %d\n", N);
  //printf("Ta %.2e\n", time_eval);
  //printf("Td %.2e\n", time_chck);
  //printf("Rt %.2e\n", time_chck/time_eval);
  //printf("Ea %.2e\n", relerr);
  printf("Run time for bfio_eval: %f\n", end_time - start_time );
  //
  return 0;
}
