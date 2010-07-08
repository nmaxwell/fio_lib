#ifndef BFIO_PROTOTYPE_CC
#define BFIO_PROTOTYPE_CC

#include "bfio_prototype.h"


int free_expansion_coefficients(double complex *** ptr, int size_x, int size_k);
double complex ** allocate_expansion_coefficients(int size_x, int size_k, int n_cheby);



int bfio(bfio_session *session)
{
    // technicalities
    
    int n_cheby = session->n_cheby;
    int N = session->N;
    int start_level = session->start_level;
    int end_level = session->end_level;
    double (*phase)(double, double, double, double) = session->phase;
    
    double unit_cheby_grid[session->n_cheby];
    for (int k=0; k<n_cheby; k++)
        unit_cheby_grid[k] = session->unit_cheby_grid[k];
    
    // initialization
    printf("Initialization.\n");
    
    int middle_level = (int)floor(0.5*(start_level+end_level));
    
    int k_level = start_level;
    int x_level = log2(N)-start_level;
    int n_kboxes = bfio_pow2(k_level);
    int n_xboxes = bfio_pow2(x_level);
    int n_boxes = bfio_pow2(start_level);
    
    double complex **next_coefficients = allocate_expansion_coefficients(n_xboxes, n_kboxes, n_cheby);
    if (next_coefficients == NULL) return __LINE__;
    
    for (int j1=0; j1<n_kboxes; j1++)
    for (int j2=0; j2<n_kboxes; j2++)
    {
        bfio_data_list *box = &(session->input_data[j1*n_boxes+j2]);
        
        if (box->size > 0)
        {
            double interp_grid1[n_cheby];
            double interp_grid2[n_cheby];
            for (int k=0; k<n_cheby; k++)
            {
                interp_grid1[k] = (unit_cheby_grid[k]+j1)/n_kboxes;
                interp_grid2[k] = (unit_cheby_grid[k]+j2)/n_kboxes;
            }
            
            for (int i1=0; i1<n_xboxes; i1++)
            for (int i2=0; i2<n_xboxes; i2++)
            {
                printf("(j1, j2), (i1, i2):\t%d\t%d\t%d\t%d\n", j1, j2, i1, i2);
                
                double x_center1 = (0.5+i1)/n_xboxes;
                double x_center2 = (0.5+i2)/n_xboxes;
                
                for (int t1=0; t1<n_cheby; t1++)
                for (int t2=0; t2<n_cheby; t2++)
                {
                    double complex temp = 0;
                    double complex value = 0;
                    for (int index=0; index < box->size; index++)
                    {
                        double source_point1 = (box->data)[index].x1;
                        double source_point2 = (box->data)[index].x2;
                        
                        temp  = lagrange_basis(source_point1, n_cheby, interp_grid1, t1);
                        temp *= lagrange_basis(source_point2, n_cheby, interp_grid2, t2);
                        temp *= cexp( +I*bfio_2pi*N*phase(x_center1, x_center2, source_point1, source_point2));
                        value += temp*((box->data)[index].value);
                    }
                    
                    value *= cexp( -I*bfio_2pi*N*phase(x_center1, x_center2, interp_grid1[t1], interp_grid2[t2]));
                    next_coefficients[n_kboxes*(n_kboxes*(n_xboxes*i1+i2)+j1)+j2][t1*n_cheby+t2] = value;
                }
            }
        }
    }
    
    double complex **coefficients = next_coefficients;
    next_coefficients = NULL;
    
    k_level -= 1;
    x_level += 1;
    
    //recurrence, first half.
    printf("Recurrence, first half.\n");
    
    while (k_level >= middle_level)
    {
        int n_kboxes = bfio_pow2(k_level);
        int n_xboxes = bfio_pow2(x_level);
        
        next_coefficients = allocate_expansion_coefficients(n_xboxes, n_kboxes, n_cheby);
        if (next_coefficients == NULL) return __LINE__;
        
        //double complex temp[n_cheby][n_cheby];
        
        for (int j1=0; j1<n_kboxes; j1++)
        for (int j2=0; j2<n_kboxes; j2++)
        for (int i1=0; i1<n_xboxes; i1++)
        for (int i2=0; i2<n_xboxes; i2++)
        {
            printf("(j1, j2), (i1, i2):\t%d\t%d\t%d\t%d\n", j1, j2, i1, i2);
            
            double interp_grid1[n_cheby];
            double interp_grid2[n_cheby];
            for (int k=0; k<n_cheby; k++)
            {
                interp_grid1[k] = (unit_cheby_grid[k]+j1)/n_kboxes;
                interp_grid2[k] = (unit_cheby_grid[k]+j2)/n_kboxes;
            }
            
            double x_center1 = (0.5+i1)/n_xboxes;
            double x_center2 = (0.5+i2)/n_xboxes;
            
            for (int t1=0; t1<n_cheby; t1++)
            for (int t2=0; t2<n_cheby; t2++)
            {
                double complex temp = 0;
                
                for (int c1=0; c1<2; c1++)
                for (int c2=0; c2<2; c2++)
                {
                    int jc1 = j1*2+c1;
                    int jc2 = j2*2+c2;
                    
                    for (int s1=0; s1<n_cheby; s1++)
                    for (int s2=0; s2<n_cheby; s2++)
                    {
                        double source_point1 = (unit_cheby_grid[s1]+jc1)/(n_kboxes*2);
                        double source_point2 = (unit_cheby_grid[s2]+jc2)/(n_kboxes*2);
                        
                        temp =  cexp( +I*bfio_2pi*N*phase(x_center1, x_center2, source_point1, source_point2));
                        temp *= lagrange_basis(source_point1, n_cheby, interp_grid1, t1);
                        temp *= lagrange_basis(source_point2, n_cheby, interp_grid2, t2);
                        temp *= coefficients[n_kboxes*(n_kboxes*(n_xboxes*(i1/2)+i2/2)+jc1)+jc2][s1*n_cheby+s2];
                    }
                }
                
                temp *= cexp( -I*bfio_2pi*N*phase(x_center1, x_center2, interp_grid1[t1], interp_grid2[t2]));
                next_coefficients[n_kboxes*(n_kboxes*(n_xboxes*i1+i2)+j1)+j2][t1*n_cheby+t2] = temp;
            }
        }
        
        free_expansion_coefficients(&coefficients, n_xboxes/2, n_kboxes*2);
        coefficients = next_coefficients;
        next_coefficients = NULL;
        
        k_level -= 1;
        x_level += 1;   
    }
    
    
    return 0;
}







int bfio_initialize_session(bfio_session *session,  int N, int n_cheby, int start_level, int end_level, double (*phase)(double, double, double, double))
{
    *session = (bfio_session){N, n_cheby, start_level, end_level, phase, NULL, 0, NULL};
    
    session->input_data_size = bfio_pow2(session->start_level*2);
    session->input_data = (bfio_data_list *)calloc(session->input_data_size, sizeof(bfio_data_list));
    if (session->input_data == NULL)
        return __LINE__;
    
    
    session->unit_cheby_grid = chebyshev_grid(n_cheby, 0.0, 1.0);
    if (session->unit_cheby_grid == NULL)
        return __LINE__;
    
    return 0;
}

int bfio_destruct_session(bfio_session *session)
{
    if (session->input_data != NULL)
    {
      int n_boxes = bfio_pow2(session->start_level);
      for (int j1=0; j1<n_boxes; j1++)
	for (int j2=0; j2<n_boxes; j2++) {
	  int error = bfio_free_input_data(session, j1, j2);
	  if (error) return __LINE__;
	}

        free(session->input_data);
        session->input_data = NULL;
    }
    
    if (session->unit_cheby_grid != NULL) {
        free(session->unit_cheby_grid);
        session->unit_cheby_grid= NULL;
      }
    
    *session = (bfio_session){0, 0, 0, 0, 0, NULL, 0, NULL};
    return 0;
}

int bfio_initialize_input_data(bfio_session *session, int j1, int j2, int size)
{
  if (session->input_data == NULL)
    return __LINE__;
    
  int n_boxes = bfio_pow2(session->start_level);
    
  if (j1>=n_boxes || j2>=n_boxes)
    return __LINE__;
    
  if (((session->input_data)[j1*n_boxes+j2]).data != NULL)
    {
      free(((session->input_data)[j1*n_boxes+j2]).data);
      ((session->input_data)[j1*n_boxes+j2]).data = NULL;
    }

  ((session->input_data)[j1*n_boxes+j2]).size = size;    
  ((session->input_data)[j1*n_boxes+j2]).data = (bfio_data_point *)calloc(size, sizeof(bfio_data_list));;
    
if (((session->input_data)[j1*n_boxes+j2]).data == NULL)
  return __LINE__;
    
return 0;
}

int bfio_free_input_data(bfio_session *session, int j1, int j2)
{
  if (session->input_data == NULL)
    return __LINE__;
  
  int n_boxes = bfio_pow2(session->start_level);

  if (j1>=n_boxes || j2>=n_boxes)
    return __LINE__;
 
  if (((session->input_data)[j1*n_boxes+j2]).data != NULL)
    {
      free(((session->input_data)[j1*n_boxes+j2]).data);
      ((session->input_data)[j1*n_boxes+j2]).data = NULL;
    }
  ((session->input_data)[j1*n_boxes+j2]).size = 0;

  return 0;
    
}

int bfio_set_input_data_array(bfio_session *session, int j1, int j2, int size, double *x1, double *x2, double complex *value)
{
  int error = bfio_initialize_input_data(session, j1, j2, size);
  if (error != 0)
    return error+10;
    
  if (x1 == NULL || x2 == NULL || value == NULL)
    return 1;
    
  int n_boxes = bfio_pow2(session->start_level);
  for (int position=0; position<size; position++)
    (((session->input_data)[j1*n_boxes+j2]).data)[position] = (bfio_data_point){x1[position], x2[position], value[position]};
    
  return 0;
}

int bfio_set_input_data_point(bfio_session *session, int j1, int j2, int position, double x1, double x2, double complex value)
{
  if (session->input_data == NULL)
    return __LINE__;
    
  int n_boxes = bfio_pow2(session->start_level);
    
  if (j1>=n_boxes || j2>=n_boxes)
    return __LINE__;
    
  if (((session->input_data)[j1*n_boxes+j2]).data == NULL)
    return __LINE__;
    
  if (position >= ((session->input_data)[j1*n_boxes+j2]).size)
    return __LINE__;
    
  (((session->input_data)[j1*n_boxes+j2]).data)[position] = (bfio_data_point){x1, x2, value};
    
    return 0;
}

int free_expansion_coefficients(double complex ***ptr, int size_x, int size_k)
{
    if (ptr != NULL)
    {
        for (int i1=0; i1<size_x; i1++)
        for (int i2=0; i2<size_x; i2++)
        for (int j1=0; j1<size_x; j1++)
        for (int j2=0; j2<size_x; j2++)
            if (ptr[size_k*(size_k*(size_x*i1+i2)+j1)+j2] != NULL)
                free(ptr[size_k*(size_k*(size_x*i1+i2)+j1)+j2]);
        
        free(ptr);
    }
    
    return 0;
}

double complex ** allocate_expansion_coefficients(int size_x, int size_k, int n_cheby)
{
    double complex **ptr = (double complex **)calloc(size_x*size_x*size_k*size_k, sizeof(double complex *));
    
    if (ptr == NULL)
        return ptr;
    
    for (int i1=0; i1<size_x; i1++)
    for (int i2=0; i2<size_x; i2++)
    for (int j1=0; j1<size_x; j1++)
    for (int j2=0; j2<size_x; j2++)
    {
        ptr[size_k*(size_k*(size_x*i1+i2)+j1)+j2] = (double complex *)calloc(n_cheby*n_cheby, sizeof(double complex));
        
        if (ptr[size_k*(size_k*(size_x*i1+i2)+j1)+j2] == NULL)
        {
            free_expansion_coefficients(&ptr, size_x, size_k);
            return NULL;
        }
    }
    
    return ptr;
}







#endif
