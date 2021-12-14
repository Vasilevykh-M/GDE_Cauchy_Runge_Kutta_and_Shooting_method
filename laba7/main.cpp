#define _USE_MATH_DEFINES
#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<cmath>

using namespace std;

const double eps1 = 0.00001; // для автоматического выбора шага
const double eps = 0.0001; // точность выполнения граничного условия

const int variant = 4;
const double ksi = 0.1;


double A(double x)
{
	return (variant == 4) ? 30 * (-x + 0.5) : 40 * (x + 1);
}

double B(double x)
{
	return (variant == 11) ? pow(x, 2) + 1 : pow(x, 2) + 2;
}

double C(double x)
{
	return (variant == 11) ? x + 2:  x + 1;
}

double Ynp(double x)
{
	return 1 + x + 10 * log(variant + 1) * pow(x, 3) * pow(1 - x, 3);
}

double DYnp_DX(double x)
{
	return 1 + 30 * log(variant + 1) * pow((x - pow(x, 2)), 2) * (1 - 2 * x);
}

double D2Ynp_DX2(double x)
{
	return 60 * log(variant + 1) * (x - pow(x, 2)) * (pow(1 - 2 * x, 2) - x + pow(x, 2));
}


double F(double x)
{
	return D2Ynp_DX2(x) + A(x) * DYnp_DX(x) - B(x) * Ynp(x) + C(x) * sin(Ynp(x));
}


double f(double x, double y, double z) /*не ебу что это*/
{
	return -A(x) * z + B(x) * y - C(x) * sin(y) + F(x);
}

double utx(double t,double x)
{
	return x + 0.1 * t * sin(M_PI * x) * variant;
}
double tauExp(double h)
{
	return pow(h, 2) / (4 * ksi);
}

double tauImp(double h)
{
	return h;
}

double f2(double t, double x)
{
	return 0.1 * sin(M_PI * x) * variant + 0.1 * 0.1 * t * variant * pow(M_PI, 2) * sin(M_PI * x);
}

void RungeKutteMethod()
{

}

vector<double> ExplicitShem(int n, int cur_i)
{
	vector<double> u_cur(n + 1);
	vector<double> u_next(n + 1);

	double h = 1. / n;
	double tau = tauExp(h);

	for (int i = 0; i < u_cur.size(); ++i)
	{
		u_cur[i] = h * i;
	}

	for (int i = 0; i < cur_i; ++i)
	{
		for (int j = 1; j < n; ++j)
		{
			u_next[j] = u_cur[j] + tau * ksi * (u_cur[j + 1] - 2 * u_cur[j] + u_cur[j - 1]) / pow(h, 2) + tau * f2(tau * i, j * h);
		}
		u_cur = u_next;
	}

	return u_cur;
}

vector<double> solveSLAU(vector<vector<double>> a, vector<double> b)
{
	vector<double> alpha(a.size());
	vector<double> beta(a.size());

	alpha[1] = -a[0][1] / a[0][0];
	beta[1] = b[0] / a[0][0];

	for (int i = 1; i < alpha.size() - 1; ++i)
	{
		alpha[i + 1] = -a[i][i + 1] / (a[i][i] + a[i][i - 1] * alpha[i]);
		beta[i + 1] = (-a[i][i - 1] * beta[i] + b[i]) / (a[i][i] + a[i][i - 1] * alpha[i]);
	}

	vector<double> x(a.size());

	x[a.size() - 1] = (b[a.size() - 1] - a[a.size() - 1][a.size() - 2] * beta[a.size() - 1]) / (a[a.size() - 1][a.size() - 1] + a[a.size() - 1][a.size() - 2] * alpha[a.size() - 1]);

	for (int i = a.size() - 2; i > -1; --i)
	{
		x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
	}

	return x;
}

vector<double> InplicitShem(int n, int cur_i)
{
	double h = 1. / n;
	double tau = tauImp(h);
	double d = tau * ksi / pow(h, 2);

	vector<vector<double>> a(n - 1,(vector<double>(n-1)));

	for (int i = 0; i < a.size(); ++i)
	{
		if (i - 1 >= 0)
			a[i][i - 1] = -d;

		a[i][i] = 1 + 2 * d;

		if (i + 1 < a.size())
			a[i][i + 1] = -d;
	}

	vector<double> u_cur(n + 1);
	for (int i = 0; i < u_cur.size(); ++i)
		u_cur[i] = i * h;

	vector<double> b(n - 1);

	for(int i=0;i<cur_i;++i)
	{
		for (int j = 0; j < b.size(); ++j)
			b[j] = u_cur[j + 1] + tau * f2(tau * (i + 1), h * (j + 1));

		b[b.size() - 1] += d;
		b = solveSLAU(a, b);

		for (int j = 0; j < b.size(); ++j)
			u_cur[j + 1] = b[j];
	}

	return u_cur;

}

void FDM()
{
	vector<double> cur_x;
	double Del_T = 0;
	for(int n=8;n<=32;n*=2)
	{
		double h = 1. / n;
		for(int i=1;i*tauExp(h)<=1;++i)
		{
			cur_x = ExplicitShem(n, i);
			double delta = 0;
			double tau = tauExp(h) * i;

			for (int j = 0; j <= n; ++j)
			{
				double max = 0;
				if (abs(utx(tau, j * h) - cur_x[i]) > max)
					max = abs(utx(tau, j * h) - cur_x[i]);
				delta = max;
			}

			if (delta > Del_T)
				Del_T = delta;
		}
	}

	for (int n = 8; n <= 32; n *= 2)
	{
		double h = 1. / n;
		for (int i = 1; i * tauImp(h) <= 1; ++i)
		{
			cur_x = InplicitShem(n, i);
			double delta = 0;
			double tau = tauImp(h) * i;

			for (int j = 0; j <= n; ++j)
			{
				double max = 0;
				if (abs(utx(tau, j * h) - cur_x[i]) > max)
					max = abs(utx(tau, j * h) - cur_x[i]);
				delta = max;
			}

			if (delta > Del_T)
				Del_T = delta;
		}
	}
}

void main()
{

}