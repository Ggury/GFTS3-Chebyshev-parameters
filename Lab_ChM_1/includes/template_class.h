#pragma once
#include <iostream>
#include <cmath>

using namespace std;

class Ax_b
{
public:
	Ax_b(double** in_Matrix, int in_n, double* in_b, double in_t_, double in_epsilon, int in_N);
	~Ax_b();

	double** Matrix_A;
	int n;                  // �� ����� ���� n*m  x(x1,x2,...,xn) y(y1,y2,...,ym)
							// n ����� ��� ��� � ����� ��� ������ ������
							// ������� �� matrix solver ������
	double* b;
	double t_;              // 1/(Mmin + Mmax) �� ������ ������. �����    (�����, ��� ����...   (TT_TT  )  )

	double epsilon;         // ��������       ����� �������� �������� ���� �� ��������� ���� �������� (��� ����� N �����)
	int N;                  // �������� ��������� - max ���-�� �����    (epsilon = 0 , ����� ����� ������� N ��������)

	//------------------------------------------------
	//     ��� ������ � ��. (����� ������, ���� �� ����� �������)

	double* r; // �������                     // r [��� ������] [��-� �������]
	double* x_0; // �����������               // x0 [��� ������] [��-� �������]

	double norma_r;
	double my_epsilon;// ������� ��������

	void S_I_M_no_input(); // x0(0,0,...,0,1)
	void S_I_M(double* x0);  // ����� ������� ��������

	void M_R_M_no_input();
	void M_R_M(double* x0); // ����� ����������� �������



private:

	
	double norma(double* x);  //  || x ||   - ����� x
	double* a_minus_b(double* a, double* b); //vectors a,b    return a-b
	double* var_x(double var, double* x); // ����� * ������
	double scalar_multiplication(double* a, double* b); //��������� ������������ (a,b)
	double* A_x(double* x); // Matrix * Vector

	

	

};

