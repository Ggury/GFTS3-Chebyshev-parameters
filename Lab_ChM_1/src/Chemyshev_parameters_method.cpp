#define M_PI 3.14159265358979323846
#include <limits>
#include "../includes/template_class.h"

void Ax_b::C_P_M_no_input()
{
	double* x0 = new double[n];

	for (int i = 0; i < n - 1; i++)
	{
		x0[i] = 0;
	}
	x0[n - 1] = 1;

	C_P_M(x0);
}


void Ax_b::C_P_M(double* x0)
{
	double min_ = std::numeric_limits<double>::max();
	double max_ = 0.;
	for (int i = 0; i < n; i++)
	{
		double sum = 0.;
		for (int j = 0; j < n; j++)
		{
			if(i!=j) sum = sum + abs(Matrix_A[i][j]);
		}
		if (sum - Matrix_A[i][i] <= min_) min_ = sum - Matrix_A[i][i];
		if (sum + Matrix_A[i][i] >= max_) max_ = sum + Matrix_A[i][i];
	}
	//cout << min_ << "    " << max_;
	double* x_prev = x0;
	r = a_minus_b(A_x(x0), b);  // r = A*x0 - b  , где x0 - начальное приближение.   Далее x0 - приближение (step)

	cout << endl << "//---[ Chebyshev Method ]---//" << endl << endl << "step=0" << endl << "  >>| ||r(0)|| = " << norma(r) << endl << endl << endl;

	int step = 1;
	//int step_for_s = 0;
	int s = 0;

	while (epsilon <= my_epsilon && step <= N)
	{
		t_ = 1. / (((min_ + max_) / 2.) + (((max_ - min_) / 2.) * cos((M_PI /	 (2. * double(k))) * (1.+ (2. * double(s))))));
		cout << endl << "t is " << t_ << endl;
		s = (step - 1) % k;

		x0 = a_minus_b(x0, var_x(t_, r));   // x0[step] = x0[step-1] - t_ * r[step-1]   , где x0 - приближение, t_ - постоянное число метода (по оценке собств. чисел), r - невязка

		r = a_minus_b(A_x(x0), b);        // r[step] = A*x0[step] - b



		norma_r = norma(r);
		my_epsilon = norma(a_minus_b(x0, x_prev));

		cout << endl << "step=" << step << endl << "  >>| ||r(" << step << ")|| = " << norma_r << endl << "  >>| epsilon(" << step << ") = " << my_epsilon << endl << " t: " << t_ << endl << endl << endl;

		x_prev = x0;
		step++;
		//step_for_s += 1;
	}
	cout << "-[ solution ]-" << endl;
	for (int i = 0; i < n; i++)
	{
		cout << "x0[" << i << "] = " << x0[i] << endl;
	}
	cout << endl << "-----------------------------------" << endl << endl;
}