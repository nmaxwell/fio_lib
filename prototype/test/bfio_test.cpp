
#include "../cppbfio/bfio_prototype.cpp"
#include <iostream>
#include <time.h>



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




int test1()
{
  int error=0;
    
  int n=1;

  for (n=0;n<=9;n++)
    for(int start_level=n; start_level>=n*3/4; start_level--)
      {
	printf("%d\t%d\n", n, start_level);
	
	int N = 1 << n;
	int n_cheby = 9;
	//int start_level = n;
	int end_level = 0;
	bfioSession session(N, n_cheby, start_level, end_level, dft_phase);
	
	srand ( time(NULL) );

	complex<double> * input = (complex<double> *)calloc(N*N, sizeof(complex<double>));
	
	for (int i=0; i<N; i++)
	  for (int j=0; j<N; j++)
	    input[i*N+j] = uniform(-1,1)+uniform(-1,1)*bfio_iu;
    
	for (int j1=0; j1<(1<<start_level); j1++)
	  for (int j2=0; j2<(1<<start_level); j2++) {
	    int size = 1<<(n-start_level);
	    for (int l1=0; l1<size; l1++)
	      for (int l2=0; l2<size; l2++) {
		double x1 = ((double)(j1*size+l1))/(1<<start_level);
		double x2 = ((double)(j2*size+l2))/(1<<start_level);
		complex<double>  value = input[(j1*size+l1)*N+j2*size+l2];
		
		session.input(j1, j2).push_back( bfioDataPoint(x1, x2, value));
	      }
	  }
	free(input);
	
	session.execute();
	
      }
  return 0;
}






int main()
{
  int error = 0;
    
  while (1)
  error = test1();
    
  if (error)
    printf("error ocurred: %d\n", error);
  else
    printf("All good.\n");
    
    
  return 0;
}
