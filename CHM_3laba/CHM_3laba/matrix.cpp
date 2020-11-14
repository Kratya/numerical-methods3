#include "matrix.h"

void SLAE::multyMatrixVector(double_v x, double_v res)
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


double SLAE::scal(double_v x, double_v y)
{
	double sum = x[0] * y[0];
	for (int i = 1; i < n; ++i)
	{
		sum += x[i] * y[i];
	}
	return sum;
}

int SLAE::LOSSolver(double eps, int maxiter)
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
	double sr = scal(r, r);
	while (sr > eps && count <= maxiter)
	{
		double pp = scal(p, p);
		double ak = scal(p, r) / pp;
		for (int i = 0; i < n; ++i)
		{
			x[i] = x[i] + ak * z[i];
			r[i] = r[i] - ak * p[i];
		}
		multyMatrixVector(r, mv);
		double bk = -scal(p, mv) / pp;
		for (int i = 0; i < n; ++i)
		{
			z[i] = r[i] + bk * z[i];
			p[i] = mv[i] + bk * p[i];
		}
		sr = sqrt(scal(r, r));
		++count;
		printf_s("%d: %.2lf", count, sr);
	}
	return count;
}

int SLAE::LOSSolverSpecifiDiag1()
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
	double normb = sqrt(scal(b, b));
	for (double sr = sqrt(scal(r, r)); sr / normb > eps && count <= maxiter; ++count)
	{
		double pp = scal(p, p);
		double ak = scal(p, r) / pp;
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
		double bk = -scal(p, mv) / pp;
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
		printf_s("%d: %.2lf", count, sr);
	}
	return count;
}



int SLAE::LOSSolverSpecifiDiag2()
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
	double normb = sqrt(scal(b, b));
	for (double sr = sqrt(scal(r, r)); sr / normb > eps && count <= maxiter; ++count)
	{
		double pp = scal(p, p);
		double ak = scal(p, r) / pp;
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
		double bk = -scal(p, mv) / pp;
		for (int i = 0; i < n; ++i)
		{
			z[i] = r[i] + bk * z[i];
			p[i] = mv[i] + bk * p[i];
		}
		sr = sqrt(scal(r, r));
		printf_s("%d: %.2lf", count, sr);
	}
	return count;
}


int SLAE::LOSSolverSpecifiDiag3()
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


	double normb = sqrt(scal(b, b));
	for (double sr = sqrt(scal(r, r)); sr / normb > eps && count <= maxiter; ++count)
	{
		double pp = scal(p, p);
		double ak = scal(p, r) / pp;
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
		double bk = -scal(p, mv) / pp;
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
		printf_s("%d: %.2lf", count, sr);
	}
	return count;
}

void SLAE::WriteSolveInFile(const char* fileName)
{
	file out;
	open_write(out, fileName);
	if (out == 0)
	{
		return;
	}
	for (int i = 0; i < n; ++i)
	{
		fprintf(out, "%.12e \n", x[i]);
	}
	fclose(out);
}



SLAE::SLAE()
{
	file in1, in2, in3, in4, in5, in6;
	bool flag = false;
	open_read(in1, "param.txt");
	open_read(in2, "ig.txt");
	open_read(in3, "jg.txt");
	open_read(in4, "gg.txt");
	open_read(in5, "di.txt");
	open_read(in6, "vector.txt");
	if (!(in1 && in2 && in3 && in4 && in5 && in6))
	{
		printf_s("Error read file!");
		throw;
	}
	fscanf_s(in1, "%d", &n);
	fscanf_s(in1, "%d", &maxiter);
	fscanf_s(in1, "%lf", &eps);
	ia = new int[n + 1];
	di = new double[n];
	b = new double[n];
	mv = new double[n];
	z = new double[n];
	r = new double[n];
	p = new double[n];
	x = new double[n];
	diag = new double[n];
	fscanf_s(in2, "%d", &ia[0]);
	if (ia[0] == 1)
	{
		flag = true;
		--ia[0];
	}
	for (int i = 1; i < n + 1; ++i)
	{
		fscanf_s(in2, "%d", &ia[i]);
		if (flag)
		{
			--ia[i];
		}
	}
	int na = ia[n];
	ja = new int[na];
	a = new double[na];
	for (int i = 0; i < na; ++i)
	{
		fscanf_s(in3, "%d", &ja[i]);
		if (flag)
		{
			--ja[i];
		}
		fscanf_s(in4, "%lf", &a[i]);
	}
	for (int i = 0; i < n; ++i)
	{
		fscanf_s(in5, "%lf", &di[i]);
	}
	for (int i = 0; i < n; ++i)
	{
		fscanf_s(in6, "%lf", &b[i]);
	}
}


SLAE::~SLAE()
{
	delete[] di;
	delete[] a;
	delete[] ia;
	delete[] ja;
	delete[] b;
	delete[] mv;
	delete[] z;
	delete[] r;
	delete[] p;
	delete[] x;
	delete[] diag;
}
