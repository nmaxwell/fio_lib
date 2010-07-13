


#include "bfio_interpolation.cpp"

#include <iostream>
#include <math.h>
#include <complex>
#include <vector>
#include <assert.h>
#include<time.h>

using namespace std;

#include <Eigen/Core>
#include <Eigen/Array>
#include <Eigen/Eigen>
USING_PART_OF_NAMESPACE_EIGEN
//#define EIGEN_NO_DEBUG

typedef Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic > MatrixXc;
typedef Matrix<complex<double>, Eigen::Dynamic, 1 > VectorXc;


template<class T >
void debug( T &X)
{

  for (int i=0; i<X.rows(); i++)
    for (int j=0; j<X.cols(); j++)
      cout << i << "\t" << j << "\t" << X(i,j) << endl;
	    
  //cout << endl << endl;
  //exit(0);
}



/*
void zgemm_(char *transa,char *transb,int *m,int *n,int *k,complex<double> *alpha,complex<double> *a,int *lda,complex<double> *b, int *ldb, complex<double> *beta, complex<double> *c, int *ldc);

int zgemm(complex<double> alpha, MatrixXc &A, MatrixXc &B, complex<double> beta, MatrixXc &C)
{
  assert( A.rows() == C.rows() );  assert(A.cols() == B.rows() );  assert( B.cols() == C.cols() );
  char transa = 'N';
  char transb = 'N';
  int m = C.rows();  int n = C.cols();  int k = A.cols();
  zgemm_(&transa, &transb, &m, &n, &k,
	 &alpha, A.data(), &m, B.data(), &k, &beta, C.data(), &m);
  return 0;
  }*/


template<class T, class S >
int lagrange_tensor(T *input, T *output, int len_input1, int len_input2, int len_output1, int len_output2, S *lagrange_matrix1, S *lagrange_matrix2)
{
  if (input == NULL || output == NULL || lagrange_matrix1 == NULL || lagrange_matrix2 == NULL) return 1;
  
  for (int i=0; i<len_output1; i++)
    for (int j=0; j<len_output2; j++)
      {
	T sum=0.0;
	for (int k=0; k<len_input1; k++)
	  for (int l=0; l<len_input2; l++)
	    sum += input[k+l*len_input1]*(lagrange_matrix1[i+k*len_output1]*lagrange_matrix2[l+j*len_input2]);
	
	output[i+j*len_output1] = sum;
      }
  
  complex<double> * temp = (complex<double> *)calloc(len_output1*len_input2, 16);

  //zgemm(1, lagrange_matrix1, input, 1, temp);
  //zgemm(1, temp, lagrange_matrix2, 1, output);
  
  free(temp);
}






class expansionCoefficients
{
public:
  int x_level;
  int k_level;
  int n_cheby;
  MatrixXc *data;

public:
  expansionCoefficients():x_level(0), k_level(0), n_cheby(0), data(NULL) {}
  expansionCoefficients(int x_level, int k_level, int n_cheby);
  ~expansionCoefficients();
  
  void operator= (expansionCoefficients &rhs);

  MatrixXc &operator() (int i1, int i2, int j1, int j2);
};

void expansionCoefficients::operator= (expansionCoefficients &rhs)
{
  if (data != NULL)
    delete [] data;
  data = NULL;
  x_level=rhs.x_level;
  k_level=rhs.k_level;
  n_cheby=rhs.n_cheby;
  data = rhs.data;
  rhs.x_level = 0;
  rhs.k_level = 0;
  rhs.n_cheby = 0;
  rhs.data = NULL;
}

MatrixXc & expansionCoefficients::operator () (int i1, int i2, int j1, int j2)
{
  int size = (1 << (2*(x_level + k_level)));
  int x_size = 1<<x_level;
  int k_size = 1<<k_level;
  int stride = n_cheby*n_cheby;
  assert(i1 < x_size and i2 < x_size);
  assert(j1 < k_size and j2 < k_size);
  assert((k_size*(k_size*(x_size*i1+i2)+j1)+j2) < size);
  assert(data != NULL);
  
  return data[(k_size*(k_size*(x_size*i1+i2)+j1)+j2)];
}

expansionCoefficients::expansionCoefficients(int x_level, int k_level, int n_cheby)
  :x_level(x_level), k_level(k_level), n_cheby(n_cheby), data(NULL) {
  int size = (1 << (2*(x_level + k_level)));
  data = new MatrixXc [size];
  assert(data != NULL);
  for (int k=0; k<size; k++)
    data[k] = MatrixXc(n_cheby, n_cheby);
}

expansionCoefficients::~expansionCoefficients()
{
  if (data != NULL)
    delete [] data;
  data = NULL;
  x_level=0;
  k_level=0;
  n_cheby=0;
}






int bfio_lexing(complex<double> *input_data, complex<double> *output_data,  int N, int start_level, int end_level, int n_cheby, double (*phase)(double, double, double ,double))
{
  if (input_data == NULL || output_data == NULL) return 1;
  
  MatrixXc input(N, N);

  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      input(i,j) = input_data[j*N+i];

  MatrixXc output = MatrixXc::Zero(N, N);

  int log2N = (int)(floor(log(N)/log(2)));
  int input_box_size = 1 << (log2N - start_level);
  int output_box_size = 1 << end_level;
  int middle_level = (start_level+end_level)/2;
  const complex<double> arg(0,6.2831853071795865);// *N

  // chebyshev grid on [0,1]

  VectorXd grid(n_cheby);
  for (int k=0; k<n_cheby; k++)
    grid[n_cheby-k-1] = (cos((double)k*3.1415926535897932384/(n_cheby-1))*0.5+0.5);


  // the following matrices are defined as the transposes of what they should be, it's backwards, but is how it was done.
  // setting up interpolation matrices from parent to child boxes
  MatrixXd tmats[2] = {MatrixXd::Zero(n_cheby, n_cheby), MatrixXd::Zero(n_cheby, n_cheby)};

  {
    VectorXd temp = grid/2;
    lagrange_matrix(tmats[0].data(), n_cheby, n_cheby, grid.data(), temp.data());
    temp = (grid/2).cwise() + 0.5;
    lagrange_matrix(tmats[1].data(), n_cheby, n_cheby, grid.data(), temp.data());
  }

  // setting up interpolation matrix from input (regularly sampled) grids, to chebyshev grids.
  MatrixXd tdir_in = MatrixXd::Zero(input_box_size, n_cheby);

  { 
    double temp[input_box_size]; 
    for (int k=0; k<input_box_size; k++)
      temp[k] = ((double)k)/input_box_size;
    lagrange_matrix(tdir_in.data(), n_cheby,  input_box_size,  grid.data(), temp);
  }

  // setting up interpolation matrix from chebyshev grids to output (regularly sampled).
  MatrixXd tdir_out = MatrixXd::Zero(output_box_size, n_cheby);

  { 
    double temp[output_box_size];
    for (int k=0; k<output_box_size; k++)
      temp[k] = ((double)k)/output_box_size;
    lagrange_matrix(tdir_out.data(), n_cheby, output_box_size, grid.data(), temp);
  }

  // transposes of these interpolation matrices
  MatrixXd mats[2] = {tmats[0].transpose(), tmats[1].transpose()};
  MatrixXd dir_in = tdir_in.transpose();
  MatrixXd dir_out = tdir_out.transpose();




  int n_zones =  1 << (log2N - start_level);
  int zone_size = N/n_zones;

  for (int z1=0; z1<n_zones; z1++)
    for (int z2=0; z2<n_zones; z2++) {
      cout << "zone: " << z1 << "\t" << z2 << endl;

      int zone_start1 = z1*zone_size;
      int zone_start2 = z2*zone_size;
      int zone_end1 = (z1+1)*zone_size;
      int zone_end2 = (z2+1)*zone_size;

      expansionCoefficients coefficients;

      for (int level=start_level; level>=middle_level; level--) {

	int k_level = level;
	int x_level = log2N - k_level;
	int n_kboxes = 1 << k_level;
	int kbox_size = N/n_kboxes;
	int n_xboxes = 1 << x_level;
	int xbox_size = N/n_xboxes;

	expansionCoefficients next_coefficients(x_level, k_level, n_cheby);

	for (int k1=zone_start1/kbox_size; k1<zone_end1/kbox_size; k1++)
	  for (int k2=zone_start2/kbox_size; k2<zone_end2/kbox_size; k2++) {

	    double k_center1 = (0.5+k1)*kbox_size;
	    double k_center2 = (0.5+k2)*kbox_size;

	    for (int x1=0; x1<n_xboxes; x1++)
	      for (int x2=0; x2<n_xboxes; x2++) {

		double x_center1 = (0.5+x1)*xbox_size;
		double x_center2 = (0.5+x2)*xbox_size;

		MatrixXc all = MatrixXc::Zero(n_cheby, n_cheby);

		if (level == start_level) {

		  // retrieve data in k-box
		  MatrixXc ext = input.block(k1*kbox_size, k2*kbox_size, kbox_size, kbox_size);

		  // apply kernel
		  for (int i=0; i<kbox_size; i++)
		    for (int j=0; j<kbox_size; j++)
		      ext(i,j) *= exp(+arg*phase(x_center1, x_center2, (k1*kbox_size+i), (k2*kbox_size+j))); // this is wrong, need to scale

		  // interpolate
		  lagrange_tensor(ext.data(), all.data(), kbox_size, kbox_size, n_cheby, n_cheby, dir_in.data(), tdir_in.data());

		}
		else {

		  for (int c1=0; c1<2; c1++)
		    for (int c2=0; c2<2; c2++)
		      {
			int kc1=2*k1+c1;
			int kc2=2*k2+c2;

			// retrieve child-coefficients
			MatrixXc ext = coefficients(x1/2, x2/2, kc1, kc2);

			// apply kernel inside
			for (int i=0; i<n_cheby; i++)
			  for (int j=0; j<n_cheby; j++)
			    ext(i,j) *= exp(+arg*phase(x_center1, x_center2, (grid[i]+kc1)*kbox_size/2, (grid[j]+kc2)*kbox_size/2)); // this is wrong, need to scale

			// interpolate
			MatrixXc temp(n_cheby, n_cheby);
			lagrange_tensor(ext.data(), temp.data(), n_cheby, n_cheby, n_cheby, n_cheby, mats[c1].data(), tmats[c2].data());
			all += temp;
		      }
		}

		// apply kernel outside
		for (int i=0; i<n_cheby; i++)
		  for (int j=0; j<n_cheby; j++)
		    all(i,j) *= exp(-arg*phase(x_center1, x_center2, (grid[i]+k1)*kbox_size, (grid[j]+k2)*kbox_size)); // this is wrong, need to scale

		//cout << "level" <<  level << "\t" << x1 << "\t" << x2 << "\t" << k1 << "\t" << k2 << endl;
		next_coefficients(x1, x2, k1, k2) = all;
	      }
	  }

	coefficients = next_coefficients;
      } //level loop, first half


      { // at middle; swithcing representations
	int level = middle_level;

	int k_level = level;
	int x_level = log2N - k_level;	
	int n_kboxes = 1 << k_level;
	int kbox_size = N/n_kboxes;
	int n_xboxes = 1 << x_level;
	int xbox_size = N/n_xboxes;

	for (int k1=zone_start1/kbox_size; k1<zone_end1/kbox_size; k1++)
	  for (int k2=zone_start2/kbox_size; k2<zone_end2/kbox_size; k2++)
	    for (int x1=0; x1<n_xboxes; x1++)
	      for (int x2=0; x2<n_xboxes; x2++) {
		//cout << x1 << "\t" << x2 <<"\t" << k1 << "\t" << k2 << endl;

		MatrixXc kernel(n_cheby*n_cheby, n_cheby*n_cheby);
		for (int i1=0; i1<n_cheby; i1++)
		  for (int i2=0; i2<n_cheby; i2++)
		    for (int j1=0; j1<n_cheby; j1++)
		      for (int j2=0; j2<n_cheby; j2++)
			kernel(i1+n_cheby*i2, j1+n_cheby*j2) = exp(+arg*phase((x1+grid[i1])*xbox_size, (x2+grid[i2])*xbox_size, (k1+grid(j1))*kbox_size, (k2+grid(j2))*kbox_size));

		coefficients(x1,x2,k1,k2).resize(n_cheby*n_cheby, 1);
		coefficients(x1,x2,k1,k2) = kernel*coefficients(x1,x2,k1,k2);
		coefficients(x1,x2,k1,k2).resize(n_cheby, n_cheby);
		//debug(coefficients(x1,x2,k1,k2));
		//cout << endl;
	      }
      }






      // recursion; second half.
      for (int level=middle_level; level>=end_level; level--) {

	int k_level = level;
	int x_level = log2N - k_level;
	int n_kboxes = 1 << k_level;
	int kbox_size = N/n_kboxes;
	int n_xboxes = 1 << x_level;
	int xbox_size = N/n_xboxes;

	expansionCoefficients next_coefficients(x_level+1, k_level-1, n_cheby);

	for (int k1=zone_start1/kbox_size; k1<zone_end1/kbox_size; k1++)
	  for (int k2=zone_start2/kbox_size; k2<zone_end2/kbox_size; k2++) {

	    double k_center1 = (0.5+k1)*kbox_size;
	    double k_center2 = (0.5+k2)*kbox_size;

	    for (int x1=0; x1<n_xboxes; x1++)
	      for (int x2=0; x2<n_xboxes; x2++) {

		//cout << level << "\t" << x1 << "\t" << x2 <<"\t" << k1 << "\t" << k2 << endl;
		//debug(coefficients(x1,x2,k1,k2));
		//cout << endl;
		double x_center1 = (0.5+x1)*xbox_size;
		double x_center2 = (0.5+x2)*xbox_size;

		MatrixXc all = coefficients(x1,x2,k1,k2);

		// apply kernel
		for (int i=0; i<n_cheby; i++)
		  for (int j=0; j<n_cheby; j++)
		    all(i,j) *= exp(-arg*phase((x1+grid[i])*xbox_size, (x2+grid[j])*xbox_size, k_center1, k_center2)); // this is wrong, need to scale

		if (level != end_level) {

		  for (int c1=0; c1<2; c1++)
		    for (int c2=0; c2<2; c2++) {

		      // interpolate
		      MatrixXc ext(n_cheby, n_cheby);
		      lagrange_tensor(all.data(), ext.data(), n_cheby, n_cheby, n_cheby, n_cheby, tmats[c1].data(), mats[c2].data());

		      // apply kernel
		      for (int i=0; i<n_cheby; i++)
			for (int j=0; j<n_cheby; j++)
			  ext(i,j) *= exp(+arg*phase((2*x1+c1+grid[i])*0.5*xbox_size, (2*x2+c2+grid[j])*0.5*xbox_size, k_center1, k_center2)); // this is wrong, need to scale
		      // sum
		      next_coefficients(2*x1+c1,2*x2+c2,k1/2,k2/2) += ext;
		    }
		}
		else {

		  // interpolate
		  MatrixXc ext(xbox_size, xbox_size);
		  lagrange_tensor(all.data(), ext.data(), n_cheby, n_cheby, xbox_size, xbox_size, tdir_out.data(), dir_out.data());
		  
		  // apply kernel
		  for (int i=0; i<xbox_size; i++)
		    for (int j=0; j<xbox_size; j++)
		      ext(i,j) *= exp(+arg*phase(i+x1*xbox_size, j+x2*xbox_size, k_center1, k_center2)); // this is wrong, need to scale

		  // debug(ext);
		  // cout << endl;
		  output.block(x1*xbox_size, x2*xbox_size, xbox_size, xbox_size) += ext;
		}


	      } // x loop
	  } // k loop

	coefficients = next_coefficients;

      } //level loop, second half

    } // zone loop

  debug(output);
  

  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      output_data[j*N+i] = output(i,j);

  return 0;
}


