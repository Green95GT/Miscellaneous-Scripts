#include<iostream>
#include <cstdio>
#include <cmath>
using namespace std;

#define	PI	3.1415926

#define is_output_time(time,i,interval)			\
	((time) >= (i)*(interval))
	
static void exact_solution(char *filename, int N, int index, double *x, double a, double time);
static void output_solution(char *filename, double *x, double *v_old, int N, int index, double time);
static double flux(double a, double dt, double dx, double vp, double v0, double vn);

static void exact_solution(char *filename, int N, int index, double *x, double a, double time)
{
	int i;
	double y;
	char fname[100];
	FILE *file;

	sprintf(fname,"%s-%d",filename,index);
	file = fopen(fname,"w");
	fprintf(file,"\"Exact soln at time = %4.2f\"\n\n",time);

	for (i = 0; i < N; ++i)                                         //Initial condictions for exact solution
	{
	    if ((x[i] - a*time) >= 1.0 && (x[i]- a*time)  <= 2.0)
		y = (sin((x[i] - a*time-1.0)*PI))*(sin((x[i] - a*time-1.0)*PI));            
	    else
		y = 0.0;
		//y=( 1.0 / sqrt( 4.0*PI*a*(time + 1.0) ))*exp( -(x[i]*x[i])/(4.0*a*(time + 1.0)) );
	    fprintf(file,"%f  %f\n",x[i],y);                            //Data file for exact solution                     
	    
	}
	fclose(file);
}

static void output_solution(char *filename, double *x, double *v_old, int N, int index, double time)
{
	char fname[100];
	FILE *file;
	int i; double v,v2;
	//v=0.0; v2=0.0;
	sprintf(fname,"%s-%d",filename,index);
	file = fopen(fname,"w");
	fprintf(file,"\"Solution at time = %4.2f\"\n\n",time);
	                                                                  //Output for scheme used
    for (i = 0; i < N; ++i)
	{
        fprintf(file,"%10.6f %10.6f\n",x[i],v_old[i]);               //Data for output of scheme used
        /*v=v+abs(v_old[i]);*/
        /*v2=v2+v*v;*/
        
    }
	fclose(file);
	/*return sqrt(v2);*/
}

static	double flux(double dx, double vp, double v0, double vn)
{
	//return (v0 - vp)/dx/2.0;
	//return (vn - v0)/dx; 
    return (vn-vp)/dx/2.0;                       
	//return 0.5*(vn + vp) - a*dt*(vn - vp)/dx/2.0;                     //Lax Friendrich Scheme for v_t+a*v_x=0
	//return (vn - 2.0*v0 + vp)/dx/dx;                                  //Used for v_t=a*v_xx
}

int main()                                       //added int p.s.      
{
	double *v_old,*v_new,*v_mid,*v_tmp,*x;
	double XL,XU;
	double a,dx,dt,time,prt_interval;
	int i,ip,N;
	//double flux();

	/* Set up mesh and allocate memory */

	N = 101;
	XL = 0.0;
	XU = 6.0;
	a = 1.0;
	dx = (XU - XL)/(N-1);
	v_old = (double*)malloc(N*sizeof(double));
	v_mid = (double*)malloc(N*sizeof(double));
	v_new = (double*)malloc(N*sizeof(double));
	x = (double*)malloc(N*sizeof(double));

	/* Initialize the solution */

	for (i = 0; i < N; ++i)
	{
	    x[i] = XL + i*dx;
	    if (x[i] >= 1.0 && x[i] <= 2.0)
		    v_old[i] = sin((x[i] - 1.0)*PI)*sin((x[i] - 1.0)*PI);
	    else
	    //v_old[i] = (1.0 / sqrt ( 4.0*PI*a )) *  exp( -(x[i]*x[i]) / (4.0*a) );
	    v_old[i]=0.0;
	    v_new[i] = v_old[i];
	    v_mid[i] = v_old[i];
	}

	time = 0.0;
	dt = 1.9*dx;
	prt_interval = 1.0;

	exact_solution("exact",N,0,x,a,time);
	output_solution("hyperb",x,v_old,N,0,time);
	ip = 1;
	while (time < 4.0)
	{
	    for (i = 1; i < N-1; ++i)
	    {
		v_new[i] = v_mid[i] - dt*a*flux(dx,v_mid[i-1],v_mid[i],v_mid[i+1]);   //Central Difference for v_t=a*v_xx
	    }                                                                         
	    time += dt;
	 
	    v_tmp = v_old;
	    v_old = v_mid;
	    v_mid = v_new;
	    v_new = v_tmp;

	    if (is_output_time(time,ip,prt_interval))
	    {
		   output_solution("hyperb",x,v_mid,N,ip,time);
           exact_solution("exact",N,ip,x,a,time);
	       ip++;
	    }
	}
}
