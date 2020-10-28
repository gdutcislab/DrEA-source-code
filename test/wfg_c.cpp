
#define _USE_MATH_DEFINES
#include <WINDOWS.H>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <malloc.h>
#include <mex.h>
#include <time.h>
using namespace std;

//#define min(x,y)  ( x < y?x:y )


void wfg(double *Output,double *Z,int M,int k, int l, int testNo, int noSols, int n);
double *s_linear(double *y, int Row, int Col ,double A);
double *b_flat(double *y, int Row, int Col, double A, double B, double C);
double *b_poly(double *y, int Row, int Col, double A);
double *r_sum(double *y, int Row, int Col, int *weights,int wei_Length);
double *h_mixed(double *x, int Row, double alpha, double A);
double *h_convex(double *x,int noSols, int uLoop);

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	int M, k, l, testNo, noSols,n;
	double *Output, *Z;
	if (nrhs == 5)
	{
		
		Z = mxGetPr(prhs[1]);
		M = (int)(*Z);
		Z = mxGetPr(prhs[2]);
		k = (int)(*Z);
		Z = mxGetPr(prhs[3]);
		l = (int)(*Z);
		Z = mxGetPr(prhs[4]);
		testNo = (int)(*Z);
		Z = mxGetPr(prhs[0]);

		noSols = mxGetM(prhs[0]);		//����
		n = mxGetN(prhs[0]);			//����


		plhs[0] = mxCreateDoubleMatrix(noSols, M, mxREAL);
		/*plhs[0] = mxCreateDoubleMatrix(noSols, M, mxREAL);*/

		Output = mxGetPr(plhs[0]);


		wfg(Output,Z, M, k, l, testNo, noSols,n);
	}
	else
	{
		mexPrintf("��Ǹ�����������������\n");
	}

}
//���󳤿����Լ���,����Z�й涨���վɣ�noSols�Ǹ�������n��ά��
void wfg(double *Output, double *Z, int M, int k, int l, int testNo, int noSols, int n)
{
	int i,ii; double *j,*jj;
	int NO_TESTS = 10;
	double D = 1;
	double *S = new double[M];
	double *F; //double *T;double *FM;
	for (i = 1; i <= M; i++)
	{
		*(S + i - 1) = (double)(2 * i);
	}
	double *A = new double[NO_TESTS * (M - 1)];	//A�Ǻ������еģ�����Y��������
	for (i = 0; i < NO_TESTS * (M - 1); i++)
	{
		*(A + i) = 1;
	}
	for (i = 1,j = (A + (M - 1) * 2); i < (M - 2); i++)
	{
		*(j + i) = 0;
	}


	double *Y = new double[noSols * n];		//���и��嶼Ҫ���˲��������⻯
	for (i = 0; i < (noSols * n); i += noSols)
	{
		double iii = (i / noSols) + 1;
		for (ii = 0; ii < noSols; ii++)
		{
			*(Y + i + ii) = (*(Z + i + ii)) / (2 * iii);
		}
	}

	//��ͬ��wfg��Ӧ��ͬ�Ĳ��Ժ���,Y����Ybar
	if (testNo == 1)
	{
		int lLoop = k + 1;
		double shiftA = 0.35;

		j = s_linear((Y + (noSols * (lLoop - 1))), noSols, (n - lLoop + 1), shiftA);
		for (i = 0,jj = (Y + (noSols * (lLoop - 1))); i < noSols * (n - lLoop + 1); i++)
		{
			*(jj + i) = *(j + i);
		}
		double biasA = 0.8; double biasB = 0.75; double biasC = 0.85;
		j = b_flat( (Y + (noSols * (lLoop - 1))), noSols, (n - lLoop + 1), biasA, biasB, biasC);

		for (i = 0, jj = (Y + (noSols * (lLoop - 1))); i < noSols * (n - lLoop + 1); i++)
		{
			*(jj + i) = *(j + i);
		}
		biasA = 0.02;
		j = b_poly(Y , noSols , n , biasA);
		for (i = 0, jj = Y; i < noSols * n; i++)
		{
			*(jj + i) = *(j + i);
		}

		//Apply fourth transformation.
		double *T = new double[noSols * M];
		//T = new double[noSols * M];
		int uLoop = M - 1;
		for (i = 0; i < uLoop; i++)
		{
			int lBnd = 1 + ((i * k) / uLoop);
			int uBnd = ((i + 1) * k) / uLoop;
			int *weights = new int[ ( uBnd - lBnd + 1) ];
			for (ii = 0; ii < (uBnd - lBnd + 1); ii++)
			{
				*(weights + ii) = (2 * (lBnd + ii));
			}
			j = r_sum((Y + (noSols * (lBnd - 1))) , noSols,(uBnd - lBnd + 1) , weights, (uBnd - lBnd + 1));

			for (ii = 0; ii < noSols; ii++)
			{
				*(T + noSols * i + ii) = *(j + ii);
			}
			delete[]weights;
		}

		int *weights = new int[(n - lLoop + 1)];
		for (ii = 0; ii < (n - lLoop + 1); ii++)
		{
			*(weights + ii) = (2 * (lLoop + ii));
		}
		j = r_sum((Y + (noSols * (lLoop - 1))), noSols, (n - lLoop + 1), weights, (n - lLoop + 1));

		for (ii = 0; ii < noSols; ii++)
		{
			*(T + noSols * (M - 1) + ii) = *(j + ii);
		}
		delete[]weights;

		//Apply degeneracy constants.
		for (i = 0; i < M - 1;i++)
		{
			for (ii = 0; ii < noSols; ii++)
			{
				*(T + i * noSols + ii) = max((*(T + i * noSols + ii)), (*(A + ((testNo - 1) * (M - 1)) + i))) * ((*(T + i * noSols + ii)) - 0.5) + 0.5;
			}
		}

		////Generate objective values.
		double *FM = h_mixed(T, noSols, 1, 5);
		F = h_convex( T , noSols , uLoop);
		for (i = 0; i < noSols;i++)
		{
			*(F + ((M - 1) * noSols) + i) = *(FM + i);
		}

		for (i = 0; i < M; i++)
		{
			for (ii = 0; ii < noSols; ii++)
			{
				*(F + i * noSols + ii) = (D * (*(T + noSols * (M - 1) + ii )) ) + ((*(S + i)) * (*(F + i * noSols + ii)));
			}
		}

	}

	//delete[]j; delete[]jj;

	for (i = 0; i < noSols * M; i++)
	{
		*(Output + i) = *(F + i);
	}

	
//for (int i = 0; i < noSols * M; i++)
//{
//	*(Output + i) = i;
//}
}

double *s_linear(double *y,int Row,int Col ,double A)
{
	int Total = Row * Col;
	double *temp = new double[Total];
	for (int i = 0; i < Total; i++)
	{
		*(temp + i) = fabs((*(y + i)) - A) / fabs( floor(A - (*(y + i)) ) + A );
	}
	return temp;
}


//Row�����У�����noSols��Col�����У�����ά��n
double *b_flat(double *y, int Row, int Col, double A, double B, double C)
{
	int Total = Row * Col;
	double * temp = new double[Total];
	//����С
	double *min1 = new double[Total];
	double *min2 = new double[Total];
	for (int i = 0; i < Total;i ++ )
	{
		*(min1 + i) = min(0 , floor((*(y + i)) - B));
		*(min2 + i) = min(0, floor( C - (*(y + i)) ));
	}

	for (int i = 0; i < Total; i++)
	{
		*(temp + i) = A + ((*(min1 + i)) * A * (B - (*(y + i))) / B) - ((*(min2 + i)) * (1 - A) * ((*(y + i)) - C) / (1 - C));
		*(temp + i) = max(0, (*(temp + i)));
	}

	delete[]min1, min2;
	return temp;
}

//double min(double a,double b)
//{
//
//}

double *b_poly(double *y, int Row, int Col, double A)
{
	int Total = Row * Col;
	double * temp = new double[Total];
	for (int i = 0; i < Total; i++)
	{
		*(temp + i) = pow((*(y + i)),A);
	}
	return temp;
}

//Row�����У�����noSols��Col�����У�����ά��n
double *r_sum(double *y, int Row, int Col, int *weights, int wei_Length)
{
	double Total = 0;
	double * temp = new double[Row * Col];
	double * ybar = new double[Row];

	for (int i = 0; i < wei_Length;i++)
	{
		Total += *(weights + i);
	}

	for (int i = 0; i < Col; i++)
	{
		for (int j = 0; j < Row; j++)
		{
			*(temp + i * Row + j) = (*(y + i * Row + j)) * (*(weights + i));
			if (i > 0)
			{
				*(temp + i * Row + j) += *(temp + (i - 1) * Row + j);
			}
		}
	}
	ybar = (temp + Row * (Col - 1));
	//delete[]temp;
	for (int i = 0; i < Row; i++)
	{
		*(ybar + i) = (*(ybar + i)) / Total;
	}

	return ybar;
}


double *h_mixed(double *x,int Row, double alpha, double A)
{
	double * f = new double[Row];
	for (int i = 0; i < Row;i++)
	{
		*(f + i) = pow((1 - *(x + i) - ( cos( (2 * A * M_PI * (*(x + i)) + (M_PI / 2))) / (2 * A * M_PI) )), alpha);
	}
	return f;
}

//noSols��x���У�uLoop��x����
double* h_convex(double *x, int noSols, int uLoop)
{
	int i, ii,iii; 
	int M = uLoop + 1;
	double *F = new double[noSols * M];
	for (i = 0; i < noSols * M; i++)
	{
		*(F + i) = 1;
	}
	//double *temp = new double[noSols];
	for (i = 0; i < noSols; i++)
	{
		
		for (ii = 0; ii < uLoop; ii++)
		{
			*(F + i) *= (1 - cos((*(x + i + ii * noSols)) * M_PI / 2));
		}
	}

	for (i = 1; i < uLoop; i++)	//��ǰ�ڼ�������ֵ
	{
		for (ii = 0; ii < noSols; ii++)//��ǰ�ڼ�������
		{
			for (iii = 0; iii < (M - i - 1); iii++)	//��ǰ�ĵڼ���xֵ(��)
			{
				*(F + i * noSols + ii) *= (1 - cos((*(x + ii + iii * noSols)) * M_PI / 2)) ;
			}
			*(F + i * noSols + ii) *= (1 - sin((*(x + ii + noSols * (M - i - 1))) * M_PI / 2));
		}
	}

	for (i = 0; i < noSols; i++)
	{
		*(F + i + (M - 1) * noSols) = 1 - sin((*(x + i)) * M_PI / 2);
	}

	return F;
	//for (i = 0; i < noSols* M;i++)
	//{
	//	*(F + i) = i;
	//}
}