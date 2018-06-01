#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
using namespace std;

#define	PI	3.1415926

void truncation_error(double *v_mid, double *x, double a, double time, int N)
{
     double q = 0.0, k = 0.0, w = 0.0;
     double *p,*g;
     int i;
     
     g = (double *)malloc(N*sizeof(double));
     p = (double *)malloc(N*sizeof(double));
     if((g || p) == 0.0)
	{
		fprintf(stderr, "out of memory\n");
	}
	         
    for (i = 0; i < N; ++i)            //Initial condictions for exact solution
    {
        if ((x[i] - a*time) >= 1.0 && (x[i]- a*time) <= 2.0)                   
            g[i] = 2.0;           
	    else
	        g[i] = 1.0;
    }
    
    for (i = 0; i < N; ++i)
    {
        p[i] = abs(v_mid[i] - g[i]);
        q += p[i]*p[i];                                                  
        k += p[i];                                                                                                                                    
        if (i >= 1)                                                                  
        {
           if (w < p[i])                                                               
              w = p[i - 1];                                                           
        }      
              
    }
    //cout<<"The Frobenius Norm is = "<<sqrt(q)<<endl;                                                            
    cout<<"The 2-norm is = "<<sqrt(q)<<endl;                                  
    cout<<"The 1-norm is = "<<k<<endl;
    cout<<"The inf-norm is = "<<w<<endl<<endl;
}     
