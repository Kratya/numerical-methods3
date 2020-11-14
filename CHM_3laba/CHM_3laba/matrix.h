#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
typedef double* double_v;
typedef int* int_v;
typedef FILE* file;
using namespace std;
#define open_read(in, str) fopen_s(&in, str, "r") 
#define open_write(out, str) fopen_s(&out, str, "w") 


class SLAE
{
private:
	double eps;
	double_v di, a;
	int_v ia, ja;
	double_v b;
	double_v x;
	double_v mv;
	double_v z;
	double_v r;
	double_v p;
	double_v diag;
	int n, maxiter;
public:
	void multyMatrixVector(double_v x, double_v res);
	double scal(double_v x, double_v y);
	int LOSSolver(double eps = 1e-12, int maxiter = 10000);
	int LOSSolverSpecifiDiag1();
	int LOSSolverSpecifiDiag2();
	int LOSSolverSpecifiDiag3();
	void WriteSolveInFile(const char* fileName);
	SLAE();
	~SLAE();
};
