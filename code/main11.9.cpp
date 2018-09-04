/*!
 * \file main.cpp
 *
 * \author Ran Luo
 * \modified by Xinyi and Si
 * \date October 2016
 *
 * All matrix is stored in one dimension array by column major.
 * 
 */


#include <iostream>
#include <fstream>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#include <cmath>





#include <iostream>
#include <fstream>


#define r 2
#define c 2

#define matrixsize 4 // rXc

using namespace std;


//GMS
void minuesVector(double  v1[r], double v2[r]);
void normlizeVector( double v[r]);
void projectVector(double v[r], double u[r], double results[r]);
double vectorInnerProduct(double v1[r], double v2[r]);
void doGramSchmidt(double input_matrix[matrixsize], double output_Matrix[matrixsize], int ifnormlize);

//LU
void doLU(double input_Matrix[matrixsize], double L[matrixsize], double U[matrixsize]);
int maxAbsValueIndex(int begin, int end, double vector[r]);
void swapRows(int row1, int row2, double inputMatrix[matrixsize]);

void minuesVector(double v1[r], double v2[r])
{
	for (int i = 0; i < r; i++)
	{
		v1[i] -= v2[i];
	}
}

void normlizeVector(double v[r])
{
	double norm = 0;
	for (int i = 0; i < r; i++)
	{
		norm += v[i] * v[i];
	}
	norm = std::sqrt(norm);

	for (int i = 0; i < r; i++)
	{
		v[i] /= norm;
	}
}


void projectVector(double v[r], double u[r], double results[r])
{
	double vu = vectorInnerProduct(v, u);
	double uu = vectorInnerProduct(u, u);
	double proj = vu / uu;

	for (int i = 0; i < r; i++)
	{
		results[i] = proj*u[i];
	}
}

double vectorInnerProduct(double v1[r], double v2[r])
{
	double result = 0;
	for (int i = 0; i < r; i++)
	{
		result += v1[i] * v2[i];
	}
	return result;
}


void doGramSchmidt(double input_matrix[matrixsize], double output_Matrix[matrixsize], int ifnormlize)
{
	double buffer[r];

	//ui
	for (int i = 0; i < c; i++)
	{
		double *ui = &output_Matrix[i*r];
		double *vi = &input_matrix[i*r];
		//memcpy(ui, vi, r*sizeof(double));

		for(int j=0;j<r;j++)
		{
			ui[j] = vi[j];
		}


		for (int j = 0; j < i; j++)
		{
			double* uj = &output_Matrix[j*r];
			projectVector(vi, uj, buffer);
			minuesVector(ui, buffer);
		}
	}
	//free(buffer);

	if (ifnormlize == 1)
	{
		for (int i = 0; i < c; i++)
		{
			double* ui = &output_Matrix[i*r];
			normlizeVector(ui);
		}
	}
}

int maxAbsValueIndex(int begin, int end, double vector[r])
{
	int maxindex = -1;
	double maxValue = vector[begin];
	for (int i = begin; i < end; i++)
	{
		if (vector[i] > maxValue)
		{
			maxValue = vector[i];
			maxindex = i;
		}
	}
	return maxindex;
}

double tempc;
int tempc_int;
#define  swapDouble(a,b) \
tempc = a; \
a = b; \
b = tempc;

#define swapInt(a,b) \
tempc_int = a; \
a = b; \
b = tempc_int;


void swapRows(int row1, int row2, double inputMatrix[matrixsize])
{
	for (int i = 0; i < c; i++)
	{
		swapDouble(inputMatrix[i*r + row1], inputMatrix[i*r + row2]);
	}
}

void doLU(double input_Matrix[matrixsize], double L[matrixsize], double U[matrixsize])
{
	int Ur, Uc, Lr, Lc;
	Ur = Lr = r;
	Uc = Lc = c;
	int row = r;
	int col = c;

	//double* bufferU = (double*)malloc(sizeof(double)*row*col);
	double bufferU[matrixsize];
	double bufferL[matrixsize];
	for(int i=0;i<matrixsize;i++)
	{
		bufferU[i] = input_Matrix[i];
		bufferL[i] = 0;
		U[i] = 0;
	}

	//memcpy(bufferU, input_Matrix, sizeof(double)*row*col);

	//int *rowindex = (int*)malloc(sizeof(int)*Lr);
	int rowindex[r];


	//double* bufferL = (double*)malloc(sizeof(double)*Lr*Lc);
	//memset(bufferL, 0, sizeof(double)*Lr*Lc);
	//memset(U, 0, sizeof(double)*Ur*Uc);

	for (int i = 0; i < Lr; i++)
	{
		rowindex[i] = i;
	}

	int minK = r;

	for(int k = 0; k < minK; k++)
	{
		double* v = &bufferU[k*row];
		int i_max = maxAbsValueIndex(k, row, v);

		if (bufferU[k*row + i_max] == 0)
		{
			//printf("matrix is singular.");
		}

		swapInt(rowindex[k], rowindex[i_max]);
		swapRows(k, i_max, bufferU);
		swapRows(k, i_max, bufferL);

		//do for all rows below pivot:
		for (int i = k + 1; i < row; i++)
		{
			double f = bufferU[k*row + i] / bufferU[k*row + k];
			for (int j = k + 1; j < col; j++)
			{
				bufferU[j*row + i] = bufferU[j*row + i] - bufferU[j*row + k] * f;
			}
			bufferU[k*row + i] = 0;
			bufferL[k*Lr + i] = f;
		}
	}


	for (int i = 0; i < minK; i++)
	{
		bufferL[i*Lr + i] = 1;
	}

	//reorder L
	for (int i = 0; i < Lr; i++)
	{
		for (int j = 0; j < Lc; j++)
		{
			L[j*Lr + rowindex[i]] = bufferL[j*Lr + i];
		}
	}

	for (int i = 0; i < Ur; i++)
	{
		for (int j = 0; j < Uc; j++)
		{
			U[j*Ur + i] = bufferU[j*row + i];
		}
	}

	//free(bufferL);
	//free(bufferU);
}


//


int main(int argc, char** argv)
{
	// this is a example

	//double inputMatrix[matrixsize];
	//double outpuMatrix[matrixsize];
	//double L[matrixsize];
	//double U[matrixsize];


	//inputMatrix[0] = 3;
	//inputMatrix[1] = 1;
	//inputMatrix[2] = 2;
	//inputMatrix[3] = 2;

	//doLU(inputMatrix, L, U);

	//for (int i = 0; i < matrixsize; i++)
	//{
	//	std::cout << U[i] << std::endl;
	//}
	
	return 0;
}

