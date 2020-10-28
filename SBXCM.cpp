
#include <WINDOWS.H>      
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <mex.h>
#include <WINDOWS.H>      
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#define PI 3.1415926535897932384626433832795029
#define EPS 1.2e-7


/* Generate a random integer. */
int irand(int);

/* Generate a random double. */
double drand(double);

void real_sbx_xover2 (double *, double *, double *, double *,double *, int);

void realmutation(double *, double *, double *, int, double);

void SBXCM(double *,  double *, double *,  double *, double *, int, double);

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) 
{
	int  dimx;
	double  *x_cld,  *xb, *x1, *x_min, *x_max, pm;

	
	dimx   = mxGetM (prhs[0]);
	xb     = mxGetPr (prhs[0]);
	x1     = mxGetPr (prhs[1]);
	x_min  = mxGetPr (prhs[2]);
	x_max  = mxGetPr (prhs[3]);	
    pm     = (double)*mxGetPr (prhs[4]);

	plhs[0]= mxCreateDoubleMatrix(dimx, 1, mxREAL);
	x_cld  = mxGetPr (plhs[0]);

	SBXCM(x_cld, xb, x1,x_min, x_max, dimx, pm);
}



void SBXCM(double *x_cld, double *xb, double *x1, double *x_min,double *x_max, int dimx, double pm)
{
	real_sbx_xover2 (x_cld, xb, x1, x_min, x_max, dimx);

    realmutation(x_cld, x_min, x_max, dimx, pm);

}


void real_sbx_xover2 (double *x_cld, double *xb, double *x1, double *x_min,double *x_max, int dimx)
{
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;
	double eta_c = 30;
    if (drand(1.0) <= 1.0) 
    {
        for (int i=0; i<dimx; i++)
        {
            if (drand(1.0)<=0.5 )
            {
                if (fabs(xb[i]-x1[i]) > EPS)
                {
                    if (xb[i] < x1[i]){
                        y1 = xb[i];
                        y2 = x1[i];
                    }
                    else{
                        y1 = x1[i];
                        y2 = xb[i];
                    }
                    yl = x_min[i];
                    yu = x_max[i];
                    rand = drand(1.0);
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(eta_c+1.0));
                    if (rand <= (1.0/alpha)){
                        betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                    }
                    else{
                        betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c1 = 0.5*((y1+y2)-betaq*(y2-y1));

                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                    alpha = 2.0 - pow(beta,-(eta_c+1.0));
                    if (rand <= (1.0/alpha)){
                        betaq = pow ((rand*alpha),(1.0/(eta_c+1.0)));
                    }
                    else{
                        betaq = pow ((1.0/(2.0 - rand*alpha)),(1.0/(eta_c+1.0)));
                    }
                    c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                    if (c1<yl)
                        c1=yl;
                    if (c2<yl)
                        c2=yl;
                    if (c1>yu)
                        c1=yu;
                    if (c2>yu)
                        c2=yu;
                    if (drand(1.0)<=0.5){
                        x_cld[i] = c2;
                    }
                    else{
                        x_cld[i] = c1;
                    }
                }
                else
                {
                    x_cld[i] = xb[i];
                }
            }
            else{
                x_cld[i] = xb[i];
            }
        }
    }
    else
    {
        for (int i=0; i<dimx; i++)
        {
            x_cld[i] = xb[i];
        }
    }
    return;
}


void realmutation(double *x_cld, double *x_min, double *x_max, int dimx, double pm)
{
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;
	double eta_m = 20;

    for (int j=0; j<dimx; j++)
    {
        if (drand(1.0)<=pm)
        {
            y  = x_cld[j];
            yl = x_min[j];
            yu = x_max[j];
            delta1 = (y-yl)/(yu-yl);
            delta2 = (yu-y)/(yu-yl);
            rnd = drand(1.0);
            mut_pow = 1.0/(eta_m+1.0);
            if (rnd <= 0.5)
            {
                xy = 1.0-delta1;
                val = 2.0*rnd+(1.0-2.0*rnd)*(pow(xy,(eta_m+1.0)));
                deltaq =  pow(val,mut_pow) - 1.0;
            }
            else
            {
                xy = 1.0-delta2;
                val = 2.0*(1.0-rnd)+2.0*(rnd-0.5)*(pow(xy,(eta_m+1.0)));
                deltaq = 1.0 - (pow(val,mut_pow));
            }
            y = y + deltaq*(yu-yl);
            if (y<yl)
                y = yl;
            if (y>yu)
                y = yu;
            x_cld[j] = y;
        }
    }
    return;
}




/* Generate a random integer. */
int irand(int range)
{
     int j;
     j=(int) ((double) range * (double) rand() / (RAND_MAX + 1.0));
     return (j);
}


/* Generate a random double. */
double drand(double range)
{
     double j;
     j=(range * (double) rand() / (RAND_MAX + 1.0));
     return (j);
}