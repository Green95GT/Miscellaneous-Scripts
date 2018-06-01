#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

using namespace std;

#define	PI	3.1415926

#define is_output_time(time,i,interval)			\
	((time) >= (i)*(interval))
	
void output_solution_sys(char *, double *, double *h, double **, int, int, double);
double flux_sys(double, double, int, double, double, double);

int main(int argc, char *argv[])
{
    double **v_old,**v_new,**v_mid,**v_tmp,*x,;
    double XL,XU;
    double dx,dh,dt,time,prt_interval;
	int i,j,ip,N;	
	
	N = 121;
	XL = 0.0;
    XU = 7.0;
    dx = (XU - XL)/(N-1);
    x = (double *)malloc(N*sizeof(double));
    v_old = (double **)malloc(N*sizeof(double *));                              //Memory allocation for matrix           
	v_mid = (double **)malloc(N*sizeof(double *));                                        
	v_new = (double **)malloc(N*sizeof(double *));
	if ((v_old || v_mid || v_new || x) == 0.0)                                                               
	{                                                                                                 
		fprintf(stderr, "out of memory\n");
	}
	for (i = 0; i < N; i++)                                                      //Still memory allocation for matrix                               
	{                                                                                                            
		v_old[i] = (double *)malloc(N*sizeof(double));                                                                         
		v_mid[i] = (double *)malloc(N*sizeof(double));                                           
		v_new[i] = (double *)malloc(N*sizeof(double));
        if ((v_old[i] || v_mid[i] || v_new[i]) == 0.0)
	   {
	       fprintf(stderr, "out of memory\n");
	   }
    }
    for (j = 0; j < 3; ++j)
	{
        for (i = 0; i < N; ++i)                                                                                                // For HW #3 2-D
        {
            x[i] = XL + i*dx;
            if (x[i] >= 3.0 && x[i] <= 4.0)                                                                                     //For HW #1,2
            {
               if (j == 0)
	              v_old[i][j] = -0.0975*sin(PI*(x[i] - 3.0)) + 0.8042*sin(2.0*PI*(x[i] - 3.0));
	              //v_old[i][j] = sin(PI*(x[i] - 3.0));
	           else if (j == 1)
                  v_old[i][j] = 1.2420*sin(PI*(x[i] - 3.0)) + 0.2318*sin(2.0*PI*(x[i] - 3.0));   
                  //v_old[i][j] = -sin(PI*(x[i] - 3.0));
	           else  
	              v_old[i][j] = 0.6692*sin(2.0*PI*(x[i] - 3.0)) - 0.5474*sin(2.0*PI*(x[i] - 3.0));
                  //v_old[i][j] = sin(2.0*PI*(x[i] - 3.0));
             }
             else
                  v_old[i][j] = 0.0;
        v_new[i][j] = v_old[i][j];                                                                                                      //??
	    v_mid[i][j] = v_old[i][j];                                                                                                      //??
        }
	}
    
    time = 0.0;
	dt = 0.625*dx;
	prt_interval = 1;
    
    output_solution_sys("hyperb",x,h,v_old,N,0,time);                                                                                   //Number goes with file index - Ex: hyperb-0, hyperb-1, ...
	ip = 1;
    
    while (time < 1.0)
	{
       for (j = 0; j < 3; ++j)                                                                                                         //For HW #3
       {
           for (i = 1; i < N-1; ++i)
	       {
               v_new[i][j] = v_mid[i][j] - flux_sys(dx,dt,j,v_mid[i-1][j],v_mid[i][j],v_mid[i+1][j]);
           }
       }
       time += dt;                                                                                                                     //advances time by dt                
	   v_tmp = v_old;                                                                                                                    //Advances solution on mesh points                                                                                              
	   v_old = v_mid;
	   v_mid = v_new;
	   v_new = v_tmp;                          
	   
	   if (is_output_time(time,ip,prt_interval))
	        {
               output_solution_sys("hyperb",x,h,v_mid,N,ip,time);
               ip++;
            }
    }
    system("PAUSE");
    return EXIT_SUCCESS;
}
