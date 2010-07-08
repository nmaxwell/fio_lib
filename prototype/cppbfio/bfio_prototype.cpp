#ifndef BFIO_PROTOTYPE_CPP
#define BFIO_PROTOTYPE_CPP

#include "bfio_prototype.h"



class expansionCoefficients
{
public:
  int x_level;
  int k_level;
  int n_cheby;
  complex<double> *data;

public:
  expansionCoefficients():x_level(0), k_level(0), n_cheby(0), data(NULL) {}
  expansionCoefficients(int x_level, int k_level, int n_cheby);
  ~expansionCoefficients();

  complex<double> * retrieve(int i1, int i2, int j1, int j2);
};

complex<double> * expansionCoefficients::retrieve(int i1, int i2, int j1, int j2)
{
  int x_size = 1<<x_level;
  int k_size = 1<<k_level;
  int stride = n_cheby*n_cheby;
  assert(i1 < x_size and i2 < x_size);
  assert(j1 < k_size and j2 < k_size);
  assert(data != NULL);
  
  return &(data[(k_size*(k_size*(x_size*i1+i2)+j1)+j2)*stride]);
}

expansionCoefficients::expansionCoefficients(int x_level, int k_level, int n_cheby)
:x_level(x_level), k_level(k_level), n_cheby(n_cheby), data(NULL) {
  int size = (1 << 2*x_level)*(1<< 2*k_level)*(n_cheby*n_cheby);
  data = (complex<double> *)calloc(size, 16);
  assert(data != NULL);
}

expansionCoefficients::~expansionCoefficients()
{
  if (data != NULL)
    free(data);
  data = NULL;
  x_level=0;
  k_level=0;
  n_cheby=0;
}







bfioSession::bfioSession():N(0), n_cheby(0), start_level(0), end_level(0), phase(NULL) {}

bfioSession::bfioSession( const bfioSession &rhs):N(rhs.N), n_cheby(rhs.N), start_level(rhs.start_level), end_level(rhs.end_level), phase(rhs.phase), unit_cheby_grid(NULL), input_data(NULL)
{
  this->initialize();
}

bfioSession::bfioSession( int N, int n_cheby, int start_level, int end_level, double(*phase)(double, double ,double, double))
  :N(N), n_cheby(n_cheby), start_level(start_level), end_level(end_level), phase(phase), unit_cheby_grid(NULL), input_data(NULL)
{
  this -> initialize();
}

int bfioSession::bfioSession_free()
{
  if (input_data != NULL)
    delete [] input_data;
  input_data = NULL;
  input_data_size = 0;
  
   if (unit_cheby_grid != NULL)
     free(unit_cheby_grid);
   unit_cheby_grid = NULL;

   return 0;
 }

 bfioSession::~bfioSession()
 {
   this->bfioSession_free();

  n_cheby = 0;
  N = 0;
  phase = NULL;
  start_level = 0;
  end_level = 0; 
}

int bfioSession::initialize()
{
  this->bfioSession_free();

  input_data_size = bfio_pow2( start_level);
  input_data = new vector<bfioDataPoint> [input_data_size*input_data_size];
  
  unit_cheby_grid = chebyshev_grid(n_cheby, 0.0, 1.0);
  
  return 0;
}

vector<bfioDataPoint>& bfioSession::input(int index1, int index2)
{
  assert( input_data != NULL and index1<input_data_size and index2 < input_data_size);
  return input_data[index1*input_data_size+ index2];
}


int bfioSession::execute()
{
  // initialization
  
  int middle_level = (int)floor(0.5*(start_level+end_level));
  complex<double> argument = complex<double> (0.0, bfio_2pi*N);
  
  int k_level = start_level;
  int x_level = log2(N)-start_level;
  int n_kboxes = bfio_pow2(k_level);
  int n_xboxes = bfio_pow2(x_level);
  //cout << "k_level: " << k_level << "\tx_level: " << x_level << endl;
  expansionCoefficients coefficients(x_level, k_level, n_cheby);
  
  for (int j1=0; j1<n_kboxes; j1++)
    for (int j2=0; j2<n_kboxes; j2++) {
      vector<bfioDataPoint> &box = input(j1, j2);
      if (box.size() > 0)
	{
	  double interp_grid1[n_cheby];
	  double interp_grid2[n_cheby];
	  for (int k=0; k<n_cheby; k++) {
	    interp_grid1[k] = (unit_cheby_grid[k]+j1)/n_kboxes;
	    interp_grid2[k] = (unit_cheby_grid[k]+j2)/n_kboxes;
	  }
            
	  for (int i1=0; i1<n_xboxes; i1++)
	    for (int i2=0; i2<n_xboxes; i2++) {
	      //printf("(j1, j2), (i1, i2):\t%d\t%d\t%d\t%d\n", j1, j2, i1, i2);
                
	      double x_center1 = (0.5+i1)/n_xboxes;
	      double x_center2 = (0.5+i2)/n_xboxes;
                
	      for (int t1=0; t1<n_cheby; t1++)
		for (int t2=0; t2<n_cheby; t2++) {
		  complex<double> temp = 0;
		  complex<double> value = 0;
		  
		  for (int index=0; index<box.size(); index++) {
		    double source_point1 = box[index].x1;
		    double source_point2 = box[index].x2;
		    
		    temp  = lagrange_basis(source_point1, n_cheby, interp_grid1, t1);
		    temp *= lagrange_basis(source_point2, n_cheby, interp_grid2, t2);
		    
		    temp *= exp( +argument*phase(x_center1, x_center2, source_point1, source_point2));
		    value += temp*(box[index].value);
		  }
		  
		  value *= exp( -argument*phase(x_center1, x_center2, interp_grid1[t1], interp_grid2[t2]));
		  (coefficients.retrieve(i1,i2,j1,j2))[t1*n_cheby+t2] = value;
		}
	    }
	}
    }
  
  bfio_core(coefficients.data, N, start_level, end_level, n_cheb, phase);

  
  
}





int bfio_core(complex<double> *coefficients, int N, int start_level, int end_level, int n_cheb, double (*phase)(double, double, double ,double))
{
  
  


  return 0;
}












#endif
