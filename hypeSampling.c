
#include <WINDOWS.H>      
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <mex.h>

double hypeSampling(double* points, int popsize, double lowerbound,
		double upperbound, int nrOfSamples, int param_k );
		
int weaklyDominates( double *point1, double *point2, int no_objectives );
double drand( double from, double to );

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) 
{
	int  i, j, popsize, param_k, nrOfSamples;
	double  *f, *x, lowerbound, upperbound;
    popsize     = mxGetN(prhs[0]);
	param_k     = popsize;
	nrOfSamples = 10000000;
	plhs[0] = mxCreateDoubleMatrix (1, popsize, mxREAL);
	double *points
	point   =mxGetPr(prhs[0]);
		for (i=0;i<popsize;i++)
	{	lowerbound(i)=point(1,i);
	     for(j=1,j<dim,j++)
			 if(point(j,i)<lowerbound(i))
				 lowerbound(i)=point(j,i);
    }
	upperbound = mxGetPr(prhs[1]);
	f=hypeSampling(double* points, int popsize, double lowerbound,
		double upperbound, int nrOfSamples, int param_k );
}


double hypeSampling(double* points, int popsize, double lowerbound,
		double upperbound, int nrOfSamples, int param_k )
/**
 * Sampling the hypeIndicator
 * \f[ \sum_{i=1}^k \left( \prod_{j=1}^{i-1} \frac{k-j}{|P|-j} \right) \frac{ Leb( H_i(a) ) }{ i } \f]
 *
 * @param[out] val vector of all indicators
 * @param[in] popsize size of the population \f$ |P| \f$
 * @param[in] lowerbound scalar denoting the lower vertex of the sampling box
 * @param[in] upperbound scalar denoting the upper vertex of the sampling box
 * @param[in] nrOfSamples the total number of samples
 * @param[in] param_k the variable \f$ k \f$
 * @param[in] points matrix of all objective values dim*popsize entries
 * @param[in] rho weight coefficients
 * @pre popsize >= 0 && lowerbound <= upperbound && param_k >= 1 &&
 * 		param_k <= popsize
 */
{
	assert(popsize >= 0 );
	assert( lowerbound <= upperbound );
	assert( param_k >= 1 );
	assert( param_k <= popsize );
	int i,j;
	double rho[param_k+1];
	/** Set alpha */
	rho[0] = 0;
	for( i = 1; i <= param_k; i++ )
	{
		rho[i] = 1.0 / (double)i;
		for( j = 1; j <= i-1; j++ )
			rho[i] *= (double)(param_k - j ) / (double)( popsize - j );
	}

	int i, s, k;
	int hitstat[ popsize ];
	int domCount;
    double* val;
	double sample[nrOfSamples][dim ];	
	for( s = 0; s < nrOfSamples; s++ )
	{
		for( k = 0; k < dim; k++ )
			sample[s][k] = drand( lowerbound, upperbound );
    }
		domCount = 0;
		for( i = 0; i < popsize; i++ )
		{
			for( k = 0; k < dim; k++ )
			if( weaklyDominates( points[i], sample[j], dim) )						
			    val[i] += rho[j];		     
		}
		for( i = 0; i < popsize; i++ )
	{
		val[i] = val[i] * pow( (upperbound[i]-lowerbound[i]), dim ) / (double)nrOfSamples;
	}
	return(val)
}

int weaklyDominates( double *point1, double *point2, int no_objectives )
{
	int better;
	int i = 0;
	better = 1;


	while( i < no_objectives && better )
	{
		better = point1[i] <= point2[i];
		i++;
	}
	return better;
}

double drand( double from, double to )
{
	double j;
	j = from + (double)( (to-from)*rand() / ( RAND_MAX + 1.0) );
	return (j);
}