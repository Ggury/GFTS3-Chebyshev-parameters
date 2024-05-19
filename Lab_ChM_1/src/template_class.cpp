#include "../includes/template_class.h"



Ax_b::Ax_b(double** in_Matrix, int in_n, double* in_b, double in_t_, double in_epsilon, int in_N, int in_k)
{
	Matrix_A = in_Matrix;
	n = in_n;

	b = in_b;
	t_ = in_t_;
	k = in_k;

	epsilon = in_epsilon;
	N = in_N;

	//------------------------------

	

	r = new double[n];
	

	my_epsilon = in_epsilon + 1.;
}

Ax_b::~Ax_b()
{
}



double Ax_b::norma(double* x)
{
	double res = 0.;
	for (int i = 0; i < n; i++)
	{
		if (abs(x[i]) > res) res = abs(x[i]);
	}
	return res;
}

void Ax_b::findM(int n, int m, double h, double k)
{
	double Mmin = std::numeric_limits<double>::max();
	double Mmax = 0;
	cout << h << " : " << k << endl;
	double h2 = -1 / (h * h);
	double k2 = -1 / (k * k); 
	cout << h2 << " : " << k2 << endl;
	double a2 = -2 * (h2 + k2);
	cout << "a: " << a2 << endl;
	n = n - 1;
	m = m - 1;
	int size = n * m;
	//double sum = 0.;
	for (int i = 0; i < size; i++)
	{
		double sum = 0.;
		sum = i - 1 >= 0 ? sum + abs(h2) : sum;
		sum = i + 1 < size ? sum + abs(h2) : sum;
		sum = i - 4 >= 0 ? sum + abs(k2) : sum;
		sum = i + 4 < size ? sum + abs(k2) : sum;
		Mmin = abs(a2) - sum < Mmin ? Mmin = abs(a2) - sum : Mmin;
		Mmax = abs(a2) + sum > Mmax ? Mmax = abs(a2) + sum : Mmax;
		cout << sum << " : ";
	}
	cout << endl <<"Mmin: " << Mmin << " Mmax "<< Mmax << endl;
}



double* Ax_b::A_x(double* x)          // Matrix * Vector
{
	double* res = new double[n];

	for (int i = 0; i < n; i++)
	{
		res[i] = 0;
		for (int j = 0; j < n; j++)
		{
			res[i] += Matrix_A[i][j] * x[j];
		}
	}
	return res;
}

double* Ax_b::a_minus_b(double* a,double* d)          // Vector - Vector
{
	double* res = new double[n];

	for (int i = 0; i < n; i++)
	{
		res[i] = a[i]-d[i];
		
	}
	return res;
}

double* Ax_b::var_x(double var,double* x)
{
	double* res = new double[n];

	for (int i = 0; i < n; i++)
	{
		res[i] = var * x[i];

	}
	return res;
}

double Ax_b::scalar_multiplication(double* a, double* b)
{
	double res = 0.;
	for (int i = 0; i < n; i++)
	{
		res += a[i] * b[i];
	}
	return res;
}



