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

//切记切记，下标数字跟matlab的不一样！！少1
//切记切记，所有的矩阵都被用来排成了一列！按照跟matlab一样的排序法，一列接一列地排！！
//所有，做操作的时候，一定要：隔popsize个数的数都是同一列的，
//例如：第一行的元素index为：0,popsize,2*popsize,3*popsize....(如果一共有popsize行)

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
	int Dim,popsize, RefSize;
	double *Output, *points, *bounds;//points应该跟matlab里面的一样，行表示个体，行数即是种群个数，列表示对应维度
							//bounds也是行向量
	
	if (nrhs == 2)
	{
		points = mxGetPr(prhs[0]);
		bounds = mxGetPr(prhs[1]);
		popsize = mxGetM(prhs[0]);		//行数
		Dim = mxGetN(prhs[0]);
		RefSize = mxGetN(prhs[1]);


		//plhs[0] = mxCreateDoubleMatrix(20, Dim, mxREAL);
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		//plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);

		Output = mxGetPr(plhs[0]);
		
		if (Dim != RefSize)
		{
			mexPrintf("抱歉，输入参数有误，BOUNDS点维数应该跟个体维数一样:个体维数：");
			//mexPrintf(mxArrayToString((const mxArray *)Dim));
			//mexPrintf(" BOUNDS点维数:"）; 
			//mexPrintf(mxArrayToString((const mxArray *)RefSize));
			//mexPrintf("\n"）;
			return;
		}

		hypeSampling(Output,points, bounds, Dim, popsize);
	}
	else
	{
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		mexPrintf("抱歉，输入参数有误\n");
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
	j = from + (double)((to - from)*rand() / (RAND_MAX + 1.0));		//这里产生的是一个在to和from之间的数
	return (j);
}

//求每一列的最大值,注意，使用该函数时，一定要对应其行列数，因为输入的不是矩阵参数，只是一个一维数组
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