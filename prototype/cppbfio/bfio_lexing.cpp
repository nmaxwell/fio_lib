


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

typedef Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic> MatrixXc;


template<class T, class S >
int lagrange_tensor(T *input, T *output, int n_input_points, int n_output_points, S *lagrange_matrix1, S *lagrange_matrix2)
{
  if (input == NULL || output == NULL || lagrange_matrix1 == NULL || lagrange_matrix2 == NULL) return 1;
  
  for (int i=0; i<n_output_points; i++)
    for (int j=0; j<n_output_points; j++)
      {
	T sum=0.0;
	for (int k=0; k<n_input_points; k++)
	  for (int l=0; l<n_input_points; l++)
	    sum += input[k*n_input_points+l]*(lagrange_matrix1[i*n_input_points+k]*lagrange_matrix2[l*n_input_points+j]);
	
	output[i*n_output_points+j] = sum;
      }
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



int bfio_lexing(complex<double> *input_data, complex<double> *output,  int N, int start_level, int end_level, int n_cheby, double (*phase)(double, double, double ,double))
{
  MatrixXc input(N, N);

  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      input(i,j) = input_data[i*N+j];

  int log2N = (int)(floor(log(N)/log(2)));
  int input_box_size = 1 << (log2N - start_level);
  int output_box_size = 1 << end_level;
  int middle_level = (start_level+end_level)/2;
  complex<double> arg(0,6.2831853071795865*N);
  
  // chebyshev grid on [0,1]

  VectorXd grid(n_cheby);
  for (int k=0; k<n_cheby; k++)
    grid[n_cheby-k-1] = (cos((double)k*3.1415926535897932384/(n_cheby-1))*0.5+0.5);

  // setting up interpolation matrices from parent to child boxes
  MatrixXd mats[2] = {MatrixXd::Zero(n_cheby, n_cheby), MatrixXd::Zero(n_cheby, n_cheby)};
  
  {
    VectorXd temp = grid/2;
    lagrange_matrix(mats[0].data(), n_cheby, n_cheby, grid.data(), temp.data());
    temp = (grid/2).cwise() + 0.5;
    lagrange_matrix(mats[1].data(), n_cheby, n_cheby, grid.data(), temp.data());
  }
  
  // setting up interpolation matrix from input (regularly sampled) grids, to chebyshev grids.
  MatrixXd dir_in = MatrixXd::Zero(n_cheby, input_box_size);
  
  { 
    double temp[input_box_size];
    for (int k=0; k<input_box_size; k++)
      temp[k] = ((double)k)/input_box_size;
    lagrange_matrix(dir_in.data(), input_box_size, n_cheby, temp, grid.data());
  }
  
  // setting up interpolation matrix from chebyshev grids to output (regularly sampled).
  MatrixXd dir_out = MatrixXd::Zero(output_box_size, n_cheby);

  { 
    double temp[output_box_size];
    for (int k=0; k<output_box_size; k++)
      temp[k] = ((double)k)/output_box_size;
    lagrange_matrix(dir_out.data(), n_cheby, output_box_size, grid.data(), temp);
  }


  // transposes of these interpolation matrices
  MatrixXd tmats[2] = {mats[0].transpose(), mats[1].transpose()};
  MatrixXd tdir_in = dir_in.transpose();
  MatrixXd tdir_out = dir_out.transpose();

  
  int n_zones =  1 << (log2N - start_level);
  int zone_size = N/n_zones;
  
  for (int z1=0; z1<n_zones; z1++)
    for (int z2=0; z2<n_zones; z2++) {
	
      int zone_start1 = z1*zone_size;
      int zone_start2 = z2*zone_size;
      int zone_end1 = (z1+1)*zone_size;
      int zone_end2 = (z2+1)*zone_size;
	
      expansionCoefficients coefficients;
	
      for (int level=start_level; level>=middle_level; level--) {
	  
	int n_kboxes = 1 << level;
	int kbox_size = N/n_kboxes;
	int n_xboxes = 1 << (log2N - level);
	int xbox_size = N/n_xboxes;
	
	int k_level = level;
	int x_level = log2N - k_level;
	
	expansionCoefficients next_coefficients(x_level, k_level, n_cheby);
	  
	for (int k1=zone_start1/zone_size; k1<zone_end1/zone_size; k1++)
	  for (int k2=zone_start1/zone_size; k2<zone_end1/zone_size; k2++) {
	      
	    double k_center1 = (0.5+k1)*kbox_size;
	    double k_center2 = (0.5+k2)*kbox_size;
	      
	    for (int x1=0; x1<n_xboxes; x1++)
	      for (int x2=0; x2<n_xboxes; x2++) {
		  
		//Vector2d x_center1 = Vector2d(0.5+x1, 0.5+x2)*xbox_size;
		double x_center1 = (0.5+x1)*xbox_size;
		double x_center2 = (0.5+x2)*xbox_size;
		  
		MatrixXc all(n_cheby, n_cheby);
		  
		if (level == start_level) {
		    
		  // retrieve data in k-box
		  MatrixXc ext = input.block(k1*kbox_size, k2*kbox_size, kbox_size, kbox_size);
		    
		  // apply kernel
		  for (int i=0; i<kbox_size; i++)
		    for (int j=0; j<kbox_size; j++)
		      ext(i,j) *= exp(+arg*phase(x_center1, x_center2, (k1*kbox_size+i), (k2*kbox_size+j))); // this is wrong, need to scale
		    
		  // interpolate
		  lagrange_tensor(ext.data(), all.data(), kbox_size, n_cheby, dir_in.data(), tdir_in.data());
		    
		}else {
		    
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
			    ext(i,j) *= exp(+arg*phase(x_center1, x_center2, (i+kc1)*kbox_size/2, (j+kc2)*kbox_size/2)); // this is wrong, need to scale
			
			// interpolate
			lagrange_tensor(ext.data(), all.data(), n_cheby, n_cheby, mats[c1].data(), tmats[c2].data());
		      }

		  // apply kernel outside
		  for (int i=0; i<n_cheby; i++)
		    for (int j=0; j<n_cheby; j++)
		      all(i,j) *= exp(-arg*phase(x_center1, x_center2, (i+k1)*kbox_size, (j+k2)*kbox_size)); // this is wrong, need to scale

		  next_coefficients(x1, x2, k1, k2) = all;
		}
	      }
	  }
	
	coefficients = next_coefficients;
      }
    }
  /*
    vector<Point2> so;
    for(int j=0; j<kB; j++)
    for(int i=0; i<kB; i++)
    so.push_back( Point2(i,j) );
    vector<Point2> ko;
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++)
    ko.push_back( Point2(grid(i)*kB, grid(j)*kB) );
    vector<Point2> co;
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++)
    co.push_back( Point2(grid(i)*kB/2, grid(j)*kB/2) );
    //
    for(int k1=k1stt/kB; k1<k1end/kB; k1++)
    for(int k2=k2stt/kB; k2<k2end/kB; k2++) {
    Point2 kc( (k1+0.5)*kB, (k2+0.5)*kB );
    for(int x1=0; x1<nx; x1++)
    for(int x2=0; x2<nx; x2++) {
    Point2 xc( (x1+0.5)*xB, (x2+0.5)*xB );
    //-------
		
    CpxNumMat all(NG,NG);		setvalue(all, cpx(0,0));
    if(ell==SL) {
    //get
    CpxNumMat ext(kB, kB);
    for(int j=0; j<kB; j++)
    for(int i=0; i<kB; i++) {
    ext(i,j) = f(i+k1*kB, j+k2*kB);
    }
    //scale
    vector<Point2> trg;		trg.push_back(xc);
    vector<Point2> src(so);
    for(int g=0; g<src.size(); g++)
    src[g] = src[g] + Point2(k1*kB, k2*kB);
    CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
    CpxNumMat sclaux(kB,kB,false,scl.data());
    for(int j=0; j<kB; j++)
    for(int i=0; i<kB; i++) {
    ext(i,j) = ext(i,j) * sclaux(i,j);
    }
    //transform
    CpxNumMat tmp(NG,kB);		setvalue(tmp,cpx(0,0));
    iC( zgemm(1, dir, ext, 0, tmp) );
    iC( zgemm(1, tmp, tdir, 1, all) );
    } else {
    int p1 = int(floor(x1/2));		int p2 = int(floor(x2/2));
    for(int a1=0; a1<2; a1++)
    for(int a2=0; a2<2; a2++) {
    int c1 = 2*k1+a1;		    int c2 = 2*k2+a2;
    //get
    CpxNumMat ext(PRE(c1,c2)(p1,p2));
    //scale
    vector<Point2> trg;		    trg.push_back(xc);
    vector<Point2> src(co);
    for(int g=0; g<src.size(); g++)
    src[g] = src[g] + Point2(c1*kB/2, c2*kB/2);
    CpxNumMat scl;		    iC( kernel(N, trg, src, scl) );
    CpxNumMat sclaux(NG,NG,false,scl.data());
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++)
    ext(i,j) = ext(i,j) * sclaux(i,j);
    //transform
    CpxNumMat tmp(NG,NG);		setvalue(tmp,cpx(0,0));
    iC( zgemm(1, mats(a1), ext, 0, tmp) );
    iC( zgemm(1, tmp, tmats(a2), 1, all) );
    }
    }
    //scale
    vector<Point2> trg;		trg.push_back(xc);
    vector<Point2> src(ko);
    for(int g=0; g<src.size(); g++)
    src[g] = src[g] + Point2(k1*kB,k2*kB);
    CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
    CpxNumMat sclaux(NG,NG,false,scl.data());
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++)
    all(i,j) = all(i,j) / sclaux(i,j);
    //put
    NOW(k1,k2)(x1,x2) = all;
    //cerr<<k1<<" "<<k2<<" "<<x1<<" "<<x2<<" "<<real(all(0,0))<<" "<<imag(all(0,0))<<endl;
    }//x1 x2
    }//k1k2
    PRE.resize(0,0);
    }//ell
    //----------------------------------------------------------------------
    if(1) {
    int ell = ML;
    int nk = pow2(ell);	int nx = N/nk;
    int kB = N/nk;	int xB = N/nx;
    //
    vector<Point2> ko;
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++)
    ko.push_back( Point2(grid(i)*kB, grid(j)*kB) );
    vector<Point2> xo;
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++)
    xo.push_back( Point2(grid(i)*xB, grid(j)*xB) );
    //
    for(int k1=k1stt/kB; k1<k1end/kB; k1++)
    for(int k2=k2stt/kB; k2<k2end/kB; k2++) {
    Point2 kc( (k1+0.5)*kB, (k2+0.5)*kB );
    for(int x1=0; x1<nx; x1++)
    for(int x2=0; x2<nx; x2++) {
    Point2 xc( (x1+0.5)*xB, (x2+0.5)*xB );
    //--------
    vector<Point2> src(ko);
    for(int g=0; g<src.size(); g++)
    src[g] = src[g] + Point2(k1*kB, k2*kB);
    vector<Point2> trg(xo);
    for(int g=0; g<trg.size(); g++)
    trg[g] = trg[g] + Point2(x1*xB, x2*xB);
    CpxNumMat evl(NG*NG,NG*NG);		iC( kernel(N, trg, src, evl) );
    CpxNumVec den(NG*NG,true,NOW(k1,k2)(x1,x2).data());
    CpxNumVec val(NG*NG,false,NOW(k1,k2)(x1,x2).data());
    iC( zgemv(1, evl, den, 0, val) );		//NOW(k1,k2)(x1,x2) = tmp;
    //cerr<<k1<<" "<<k2<<" "<<x1<<" "<<x2<<" "<<real(val(0))<<" "<<imag(val(0))<<endl;
    }
    }
    }//ell
    //----------------------------------------------------------------------
    for(int ell=ML; ell>=EL; ell--) {
    int nk = pow2(ell);	int nx = N/nk;
    int kB = N/nk;	int xB = N/nx;
    NumMat< NumMat< CpxNumMat > > NXT;
    if(ell!=EL) {
    NXT.resize(nk/2,nk/2);
    for(int k1=k1stt/(2*kB); k1<k1end/(2*kB); k1++)
    for(int k2=k2stt/(2*kB); k2<k2end/(2*kB); k2++) {
    NXT(k1,k2).resize(2*nx,2*nx);
    for(int x1=0; x1<2*nx; x1++)
    for(int x2=0; x2<2*nx; x2++) {
    NXT(k1,k2)(x1,x2).resize(NG,NG);		setvalue(NXT(k1,k2)(x1,x2),cpx(0,0));
    }
    }
    }
    //
    vector<Point2> to;
    for(int j=0; j<xB; j++)
    for(int i=0; i<xB; i++)
    to.push_back( Point2(i,j) );
    vector<Point2> xo;
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++)
    xo.push_back( Point2(grid(i)*xB, grid(j)*xB) );
    vector<Point2> co;
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++)
    co.push_back( Point2(grid(i)*xB/2, grid(j)*xB/2) );
    //
    for(int k1=k1stt/kB; k1<k1end/kB; k1++)
    for(int k2=k2stt/kB; k2<k2end/kB; k2++) {
    Point2 kc( (k1+0.5)*kB, (k2+0.5)*kB );
    for(int x1=0; x1<nx; x1++)
    for(int x2=0; x2<nx; x2++) {
    Point2 xc( (x1+0.5)*xB, (x2+0.5)*xB );
    //-------
    //get
    CpxNumMat all(NOW(k1,k2)(x1,x2));
    //scale
    vector<Point2> src;		src.push_back(kc);
    vector<Point2> trg(xo);
    for(int g=0; g<trg.size(); g++)
    trg[g] = trg[g] + Point2(x1*xB,x2*xB);
    CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
    CpxNumMat sclaux(NG,NG,false,scl.data());
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++) {
    all(i,j) = all(i,j) / sclaux(i,j);
    }
    //
    if(ell!=EL) {
    int q1 = int(floor(k1/2));		  int q2 = int(floor(k2/2));
    for(int a1=0; a1<2; a1++)
    for(int a2=0; a2<2; a2++) {
    int c1 = 2*x1+a1;		      int c2 = 2*x2+a2;
    //transform
    CpxNumMat ext(NG,NG);		      setvalue(ext,cpx(0,0));
    CpxNumMat tmp(NG,NG);		      setvalue(tmp,cpx(0,0));
    iC( zgemm(1, tmats(a1), all, 0, tmp) );
    iC( zgemm(1, tmp, mats(a2), 1, ext) );
    //scale
    vector<Point2> src;		      src.push_back(kc);
    vector<Point2> trg(co);
    for(int g=0; g<trg.size(); g++)
    trg[g] = trg[g] + Point2(c1*xB/2,c2*xB/2);
    CpxNumMat scl;		      iC( kernel(N, trg, src, scl) );
    CpxNumMat sclaux(NG,NG,false,scl.data());
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++)
    ext(i,j) = ext(i,j) * sclaux(i,j);
    //put
    for(int j=0; j<NG; j++)
    for(int i=0; i<NG; i++)
    NXT(q1,q2)(c1,c2)(i,j) = NXT(q1,q2)(c1,c2)(i,j) + ext(i,j);
    }
    } else {
    //transform
    CpxNumMat ext(xB,xB);		  setvalue(ext,cpx(0,0));
    CpxNumMat tmp(xB,NG);		  setvalue(tmp,cpx(0,0));
    iC( zgemm(1, tdir, all, 0, tmp) );
    iC( zgemm(1, tmp, dir, 1, ext) );
    //scale
    vector<Point2> src;		  src.push_back(kc);
    vector<Point2> trg(to);
    for(int g=0; g<trg.size(); g++)
    trg[g] = trg[g] + Point2(x1*xB, x2*xB);
    CpxNumMat scl;		  iC( kernel(N, trg, src, scl) );
    CpxNumMat sclaux(xB,xB,false,scl.data());
    for(int j=0; j<xB; j++)
    for(int i=0; i<xB; i++)
    ext(i,j) = ext(i,j) * sclaux(i,j);
    //cerr<<k1<<" "<<k2<<" "<<x1<<" "<<x2<<" "<<real(ext(0,0))<<" "<<imag(ext(0,0))<<endl;
    //put
    for(int j=0; j<xB; j++)
    for(int i=0; i<xB; i++)
    u(i+x1*xB,j+x2*xB) = u(i+x1*xB,j+x2*xB) + ext(i,j);
    }
    }//x1x2
    }//k1k2
    NOW = NXT;
    NXT.resize(0,0);
    }//ell
    }//z1z2
  


  */

  


  return 0;
}


