#include "bfio.hpp"
#include "serialize.hpp"

using namespace std;

int optionsCreate(int argc, char** argv, map<string,string>& options)
{
  options.clear();
  //2. get extra data
  for(int k=1; k<argc; k=k+2) {
    options[ string(argv[k]) ] = string(argv[k+1]);
  }
  return 0;
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
  iC( bfio.setup(opts) );
  //
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
  //
  double relerr = 0;
  int NC = 64;
  ck0 = clock();
  iC( bfio.check(f,u,NC,relerr) );
  ck1 = clock();
  double time_chck = double(ck1-ck0)/CLOCKS_PER_SEC * N * N / double(NC);
  //
  printf("RESULT\n");
  printf("N  %d\n", N);
  printf("Ta %.2e\n", time_eval);
  printf("Td %.2e\n", time_chck);
  printf("Rt %.2e\n", time_chck/time_eval);
  printf("Ea %.2e\n", relerr);
  //
  return 0;
}
