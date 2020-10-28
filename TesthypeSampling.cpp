#include <WINDOWS.H>      
#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <mex.h>
#include <WINDOWS.H>


double drand(double from, double to);
double* GetColMax(double* points, int Rows, int Cols);
double* GetColMin(double* points, int Rows, int Cols);
//double* hypeSampling(double* points, double* bounds, int Dim, int popsize);
void hypeSampling(double*Output, double* points, double* bounds, int Dim, int popsize);

//�м��мǣ��±����ָ�matlab�Ĳ�һ��������1
//�м��мǣ����еľ��󶼱������ų���һ�У����ո�matlabһ�������򷨣�һ�н�һ�е��ţ���
//���У���������ʱ��һ��Ҫ����popsize������������ͬһ�еģ�
//���磺��һ�е�Ԫ��indexΪ��0,popsize,2*popsize,3*popsize....(���һ����popsize��)

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	int Dim,popsize, RefSize;
	double *Output, *points, *bounds;//pointsӦ�ø�matlab�����һ�����б�ʾ���壬����������Ⱥ�������б�ʾ��Ӧά��
							//boundsҲ��������
	
	if (nrhs == 2)
	{
		points = mxGetPr(prhs[0]);
		bounds = mxGetPr(prhs[1]);
		popsize = mxGetM(prhs[0]);		//����
		Dim = mxGetN(prhs[0]);
		RefSize = mxGetN(prhs[1]);


		//plhs[0] = mxCreateDoubleMatrix(20, Dim, mxREAL);
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		//plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);

		Output = mxGetPr(plhs[0]);
		
		if (Dim != RefSize)
		{
			mexPrintf("��Ǹ�������������BOUNDS��ά��Ӧ�ø�����ά��һ��:����ά����");
			//mexPrintf(mxArrayToString((const mxArray *)Dim));
			//mexPrintf(" BOUNDS��ά��:"��; 
			//mexPrintf(mxArrayToString((const mxArray *)RefSize));
			//mexPrintf("\n"��;
			return;
		}

		hypeSampling(Output,points, bounds, Dim, popsize);
	}
	else
	{
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		mexPrintf("��Ǹ�������������\n");
	}
}


void hypeSampling(double* Output,double* points, double* bounds, int Dim, int popsize)
{
	int nrOfSamples = 10000000;
	//int nrOfSamples = 1000;
	double* F = mxGetPr(mxCreateDoubleMatrix(1, popsize, mxREAL));
	double* alpha = mxGetPr(mxCreateDoubleMatrix(1, popsize, mxREAL));

	for (int i = 0; i < popsize;i ++ )
	{
		alpha[i] = 1.0 / (double)(i + 1);
	}

	double *BoxL = GetColMin(points, popsize, Dim);
	double *S = (double *)malloc(nrOfSamples*Dim*sizeof(double));
	for (int i = 0; i < Dim;i++)
	{
		for (int j = 0; j < nrOfSamples; j++)
		{
			S[i * nrOfSamples + j] = drand(BoxL[i],bounds[i]);

			//Output[i * nrOfSamples + j] = S[i * nrOfSamples + j];
		}
	}
	
	double *dominated = mxGetPr(mxCreateDoubleMatrix(nrOfSamples, 1, mxREAL));

	for (int j = 0; j < popsize;j++)
	{
		for (int i = 0; i < nrOfSamples;i ++)
		{
			bool IsDominate = true;
			for (int k = 0; k < Dim; k++)
			{
				if (S[k * nrOfSamples + i] < points[k * popsize + j])
				{
					IsDominate = false;
				}
			}
			if (IsDominate == true)
			{
				dominated[i] += 1;
			}
		}
	}


	for (int j = 0; j < popsize; j++)
	{
		for (int i = 0; i < nrOfSamples; i++)
		{
			bool IsDominate = true;
			for (int k = 0; k < Dim; k++)
			{
				if (S[k * nrOfSamples + i] < points[k * popsize + j])
				{
					IsDominate = false;
				}
			}
			if (IsDominate == true)
			{
				F[j] += alpha[(int)(dominated[i] - 1)];
			}
		}
	}



	double Sequare = 1;
	for (int i = 0; i < Dim; i++)
	{
		Sequare *= (bounds[i] - BoxL[i]);
	}


	/*double *test = mxGetPr(mxCreateDoubleMatrix(popsize, 1, mxREAL));
	for (int i = 0; i < nrOfSamples;i++)
	{
		test[(int)dominated[i]] += 1;
	}*/



	for (int i = 0; i < popsize; i++)
	{
		F[i] = (F[i] * Sequare) / nrOfSamples;
		Output[0] += F[i];
	}


	//Output[0] = S[0];
	//Output = S;
}




double drand(double from, double to)
{
	double j;
	j = from + (double)((to - from)*rand() / (RAND_MAX + 1.0));		//�����������һ����to��from֮�����
	return (j);
}

//��ÿһ�е����ֵ,ע�⣬ʹ�øú���ʱ��һ��Ҫ��Ӧ������������Ϊ����Ĳ��Ǿ��������ֻ��һ��һά����
double* GetColMax(double* points, int Rows,int Cols)
{
	double *Max;
	Max = (double *)malloc(Cols*sizeof(double));
	for (int i = 0; i < Cols; i++)
	{
		for (int j = 0; j < Rows; j++)
		{
			if (j == 0)
			{
				Max[i] = points[(i * Rows) + j];
			}
			else if (Max[i] < points[(i * Rows) + j])
			{
				Max[i] = points[(i * Rows) + j];
			}
		}
	}
	return Max;
}

double* GetColMin(double* points, int Rows, int Cols)
{
	double *Min;
	Min = (double *)malloc(Cols*sizeof(double));
	for (int i = 0; i < Cols; i++)
	{
		for (int j = 0; j < Rows;j++)
		{
			if (j == 0)
			{
				Min[i] = points[(i * Rows) + j];
			}
			else if (Min[i] > points[(i * Rows) + j])
			{
				Min[i] = points[(i * Rows) + j];
			}
		}
	}
	return Min;
}