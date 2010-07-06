
#include "../cbfio/bfio_prototype.h"
#include <stdio.h>
#include <stdlib.h>
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
    
    int n=6;
    int N = 1 << n;
    int n_cheby = 9;
    int start_level = n-4;
    int end_level = 0;
    
    bfio_session session;
    error = bfio_initialize_session(&session, N, n_cheby, start_level, end_level, dft_phase);
    if (error) return error;
    
    srand ( time(NULL) );
    double complex input[N][N];
    for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
        input[i][j] = uniform(-1,1)+uniform(-1,1)*I;
    
    
    for (int j1=0; j1<(1<<start_level); j1++)
    for (int j2=0; j2<(1<<start_level); j2++)
    {
        printf("%d\t%d\n", j1, j2);
        for (int p=0; p<(session.input_data_size); p++)
            printf( "\t%d\t%d\t%lx\n", p, session.input_data[p].size, (unsigned long)(session.input_data[p].data) );
        
        
        int size = 1<<(n-start_level);
        
        
        error = bfio_intialize_input_data(&session, j1, j2, size*size);
        if (error) return error;
        
        for (int l1=0; l1<size; l1++)
        for (int l2=0; l2<size; l2++)
        {
            double x1 = ((double)(j1*size+l1))/(1<<start_level);
            double x2 = ((double)(j2*size+l2))/(1<<start_level);
            double complex value = input[j1*size+l1][j2*size+l2];
            
            int position = l1*size+l2;
            //error = bfio_set_input_data_point(&session, j1, j2, position, x1, x2, value);
            
            (((session.input_data)[j1*(1<<start_level)+j2]).data)[position] = (bfio_data_point){x1, x2, value};
            
            
            if (error) return error;
        }
        
        
    }
    
    
    
    //error = bfio(&session);
    if (error) return error;
    
    
    
    bfio_destruct_session(&session);
    if (error) return error;
    
    return 0;
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
