#pragma once
#include <stdio.h>
#include <iostream>
#include<fstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

template<class mytype>
class ALG
{
private:
	double eps;
	vector<mytype> di, a;
	vector<int> ia, ja;
	vector<mytype> b;
	vector<mytype> x;
	vector<mytype> mv;
	vector<mytype> z;
	vector<mytype> r;
	vector<mytype> p;
	vector<mytype> diag;
	int n, maxiter;
public:
	void ALG_ini();
	void multyMatrixVector(vector<mytype> &x, vector<mytype> &res);
	mytype scal(vector<mytype> &x, vector<mytype> &y);
	int LOSSolver(mytype eps = 1e-12, int maxiter = 10000);
	int LOSSolverSpecifiDiag1();
	int LOSSolverSpecifiDiag2();
	int LOSSolverSpecifiDiag3();
	void WriteSolveInFile(const char* fileName);
};

template <typename mytype>
void ALG<mytype>::multyMatrixVector(vector<mytype> &x, vector<mytype> &res)
{
	for (int i = 0; i < n; i++)
	{
		res[i] = x[i] * di[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			res[i] += a[k] * x[j];
			res[j] += a[k] * x[i];
		}
	}
}

template <typename mytype>
mytype ALG<mytype>::scal(vector<mytype> &x, vector<mytype> &y)
{
	mytype sum = x[0] * y[0];
	for (int i = 1; i < n; ++i)
	{
		sum += x[i] * y[i];
	}
	return sum;
}

template <typename mytype>
int ALG<mytype>::LOSSolver(mytype eps, int maxiter)
{
	int count = 0;
	for (int i = 0; i < n; ++i)
	{
		x[i] = 1;
	}
	multyMatrixVector(x, mv);
	for (int i = 0; i < n; ++i)
	{
		r[i] = b[i] - mv[i];
		z[i] = r[i];
	}
	multyMatrixVector(z, p);
	mytype sr = scal(r, r);
	while (sr > eps && count <= maxiter)
	{
		mytype pp = scal(p, p);
		mytype ak = scal(p, r) / pp;
		for (int i = 0; i < n; ++i)
		{
			x[i] = x[i] + ak * z[i];
			r[i] = r[i] - ak * p[i];
		}
		multyMatrixVector(r, mv);
		mytype bk = -scal(p, mv) / pp;
		for (int i = 0; i < n; ++i)
		{
			z[i] = r[i] + bk * z[i];
			p[i] = mv[i] + bk * p[i];
		}
		sr = sqrt(scal(r, r));
		++count;
		//printf_s("%d: %.2lf", count, sr);
	}
	return count;
}

template <typename mytype>
int ALG<mytype>::LOSSolverSpecifiDiag1()
{
	int count = 0;
	for (int i = 0; i < n; ++i)
	{
		x[i] = 1;
		diag[i] = sqrt(di[i]);
	}

	multyMatrixVector(x, mv);
	for (int i = 0; i < n; ++i)
	{
		r[i] = (b[i] - mv[i]) / diag[i];
		z[i] = r[i] / diag[i];
	}
	multyMatrixVector(z, p);

	for (int i = 0; i < n; ++i)
	{
		p[i] /= diag[i];
	}
	mytype normb = sqrt(scal(b, b));
	for (mytype sr = sqrt(scal(r, r)); sr / normb > eps && count <= maxiter; ++count)
	{
		mytype pp = scal(p, p);
		mytype ak = scal(p, r) / pp;
		for (int i = 0; i < n; ++i)
		{
			x[i] = x[i] + ak * z[i];
			r[i] = r[i] - ak * p[i];
		}
		for (int i = 0; i < n; ++i)
		{
			r[i] /= diag[i];
		}
		multyMatrixVector(r, mv);
		for (int i = 0; i < n; ++i)
		{
			mv[i] /= diag[i];
		}
		mytype bk = -scal(p, mv) / pp;
		for (int i = 0; i < n; ++i)
		{
			z[i] = r[i] + bk * z[i];
			p[i] = mv[i] + bk * p[i];
		}
		for (int i = 0; i < n; ++i)
		{
			r[i] *= diag[i];
		}
		sr = sqrt(scal(r, r));
		//printf_s("%d: %.2lf", count, sr);
	}
	return count;
}


template <typename mytype>
int ALG<mytype>::LOSSolverSpecifiDiag2()
{
	int count = 0;
	for (int i = 0; i < n; ++i)
	{
		x[i] = 1;
	}

	multyMatrixVector(x, mv);
	for (int i = 0; i < n; ++i)
	{
		r[i] = (b[i] - mv[i]) / di[i];
		z[i] = r[i];
	}
	multyMatrixVector(z, p);

	for (int i = 0; i < n; ++i)
	{
		p[i] /= di[i];
	}
	mytype normb = sqrt(scal(b, b));
	for (mytype sr = sqrt(scal(r, r)); sr / normb > eps && count <= maxiter; ++count)
	{
		mytype pp = scal(p, p);
		mytype ak = scal(p, r) / pp;
		for (int i = 0; i < n; ++i)
		{
			x[i] = x[i] + ak * z[i];
			r[i] = r[i] - ak * p[i];
		}
		multyMatrixVector(r, mv);
		for (int i = 0; i < n; ++i)
		{
			mv[i] /= di[i];
		}
		mytype bk = -scal(p, mv) / pp;
		for (int i = 0; i < n; ++i)
		{
			z[i] = r[i] + bk * z[i];
			p[i] = mv[i] + bk * p[i];
		}
		sr = sqrt(scal(r, r));
		//printf_s("%d: %.2lf", count, sr);
	}
	return count;
}

template <typename mytype>
int ALG<mytype>::LOSSolverSpecifiDiag3()
{
	int count = 0;
	for (int i = 0; i < n; ++i)
	{
		x[i] = 1;
	}

	multyMatrixVector(x, mv);
	for (int i = 0; i < n; ++i)
	{
		r[i] = b[i] - mv[i];
		z[i] = r[i] / di[i];
	}
	multyMatrixVector(z, p);


	mytype normb = sqrt(scal(b, b));
	for (mytype sr = sqrt(scal(r, r)); sr / normb > eps && count <= maxiter; ++count)
	{
		mytype pp = scal(p, p);
		mytype ak = scal(p, r) / pp;
		for (int i = 0; i < n; ++i)
		{
			x[i] = x[i] + ak * z[i];
			r[i] = r[i] - ak * p[i];
		}
		for (int i = 0; i < n; ++i)
		{
			r[i] /= di[i];
		}
		multyMatrixVector(r, mv);
		mytype bk = -scal(p, mv) / pp;
		for (int i = 0; i < n; ++i)
		{
			z[i] = r[i] + bk * z[i];
			p[i] = mv[i] + bk * p[i];
		}
		for (int i = 0; i < n; ++i)
		{
			r[i] *= di[i];
		}
		sr = sqrt(scal(r, r));
		//printf_s("%d: %.2lf", count, sr);
	}
	return count;
}

template <typename mytype>
void ALG<mytype>::WriteSolveInFile(const char* fileName)
{
	ofstream fout;
	fout.open(fileName);
	fout.precision(15);

	for (int i = 0; i < n; ++i)
	{
		fout << x[i] << endl;
	}
	fout.close();
}

template <typename mytype>
void ALG<mytype>::ALG_ini()
{
	ifstream in1, in2, in3, in4, in5, in6;
	bool flag = false;

	in1.open("param.txt");
	in2.open("ig.txt");
	in3.open("jg.txt");
	in4.open("gg.txt");
	in5.open("di.txt");
	in6.open("vector.txt");

	in1 >> n >> maxiter >> eps;
	
	ia.resize(n + 1);
	di.resize(n);
	b.resize(n);
	mv.resize(n);
	z.resize(n);
	r.resize(n);
	p.resize(n);
	x.resize(n);
	diag.resize(n);

	in2 >> ia[0];
	cout << "ia[0] = " << ia[0] << endl;
	if (ia[0] == 1)
	{
		flag = true;
		--ia[0];
	}
	for (int i = 1; i < n + 1; ++i)
	{
		in2 >> ia[i];
		cout << "ia[" << i << "] = " << ia[i] << endl;
		if (flag)
			--ia[i];
	}
	int na = ia[n];
	cout << "ia[n] = " << ia[n] << endl;
	ja.resize(na);
	a.resize(na);
	for (int i = 0; i < na; ++i)
	{
		in3 >> ja[i];
		cout << "ja[" << i << "] = " << ja[i] << endl;
		if (flag)
		{
			--ja[i];
		}
		in4 >> a[i];
		cout << "a[" << i << "] = " << a[i] << endl;
	}
	for (int i = 0; i < n; ++i)
	{
		in5 >> di[i];
		cout << "di[" << i << "] = " << di[i] << endl;
	}
	for (int i = 0; i < n; ++i)
	{
		in6 >> b[i];
		cout << "b[" << i << "] = " << b[i] << endl;
	}
	in1.close();
	in2.close();
	in3.close();
	in4.close();
	in5.close();
	in6.close();
}