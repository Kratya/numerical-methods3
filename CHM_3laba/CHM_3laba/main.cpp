#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include "matrix.h"

int main()
{
	ALG<double> m;
	int iter;

	clock_t s1, s2;

	m.ALG_ini();
	
	s1 = clock();
	iter = m.LOSSolver();
	s2 = clock();
	//printf_s("%d\n", s2 - s1);	
	cout << "1. Nums of iters = " <<  iter << endl;
	m.WriteSolveInFile("RESULT1.TXT");
	
	s1 = clock();
	iter = m.LOSSolverSpecifiDiag1();
	s2 = clock();
	//printf_s("%d\n", s2 - s1);
	cout << "2. Nums of iters = " << iter << endl;
	m.WriteSolveInFile("RESULT1_1.TXT");

	s1 = clock();
	iter = m.LOSSolverSpecifiDiag2();
	s2 = clock();
	//printf_s("%d\n", s2 - s1);
	cout << "3. Nums of iters = " << iter << endl;
	m.WriteSolveInFile("RESULT1_2.TXT");
/*
	s1 = clock();
	iter = m.LOSSolverSpecifiDiag3();
	s2 = clock();
	//printf_s("%d\n", s2 - s1);
	cout << "4. Nums of iters = " << iter << endl;
	m.WriteSolveInFile("RESULT1_3.TXT");

*/
	system("pause");

}
