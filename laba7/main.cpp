#define _USE_MATH_DEFINES
#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<cmath>

using namespace std;

const double eps1 = 0.00001; // дл€ автоматического выбора шага
const double eps = 0.0001; // точность выполнени€ граничного услови€

const int variant = 4;
const double ksi = 0.1;


double A(double x)
{
	return (variant == 4) ? 30 * (-x + 0.5) : 40 * (x + 1);
}

double B(double x)
{
	return (variant == 4) ? pow(x, 2) + 1 : pow(x, 2) + 2;
}

double C(double x)
{
	return (variant == 4) ? x + 2 : x + 1;
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

double H(double h, double x, double y, double& z) /*€ вроде понимаю как это работает*/
{
	double k1, k2, k3, k4;
	double l1, l2, l3, l4;

	k1 = h * z;
	l1 = h * f(x, y, z);
	k2 = h * (z + l1 / 2.);
	l2 = h * f(x + h / 2., y + k1 / 2., z + l1 / 2.);
	k3 = h * (z + l2 / 2.);
	l3 = h * f(x + h / 2., y + k2 / 2., z + l2 / 2.);
	k4 = h * (z + l3);
	l4 = h * f(x + h, y + k3, z + l3);

	y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
	z = z + (l1 + 2 * l2 + 2 * l3 + l4) / 6.;

	return y;
}


double RungeKutteMethod(double alpha, double& max_error, bool printT)
{
	double x = 0;
	double y = 1;
	double z =alpha;

	double y_prev, z_prev;
	double k1, k2, k3, k4; // дл€ уравнени€ y' = z
	double l1, l2, l3, l4; // дл€ уравнени€ z' = -Az + By - —sin(y) + F

	if (printT)
		cout << "     x     |   y(x)   |   Ypr(x)  |   z(x)   |    Delta" << endl;

	while (x < 1) { // € пон€ти€ не имею зачем нужен внешний цикл
		double h = 0.2;
		double y1, y2, y3;
		double z1;
		y_prev = y;
		z_prev = z;

		/*€ пон€л как это работает*/
		do
		{
			y_prev = y;
			z_prev = z;
			z1 = z;

			h = h / 2;

			y1 = H(h, x, y_prev, z1);
			y2 = H(h / 2., x, y_prev, z_prev);
			y3 = H(h / 2.,  x, y2, z_prev);
		} while (abs(y1 - y3) > eps1);

		double error = abs(Ynp(x) - y);
		if (error > max_error)
			max_error = error;

		if (printT)
			cout << fixed << setprecision(5) << setw(10) << x << setw(10) << y << setw(10) << Ynp(x) << setw(10) << z
			<< setw(14) << scientific << error << endl;


		y_prev = y;
		y = H(h, x, y, z);
		x += h;
	}

	return y_prev;
}


void shooting_metod()
{
	double y;
	double alpha1 = 0, alpha2 = 3, alpha;
	double B = 2;
	int itr = 0;
	double max_error;
	cout << "  itr  |   z(0)     |    y(1)     |     Delta" << endl;
	/*тут впринципе все просто*/
	do
	{
		itr++;
		alpha = (alpha2 + alpha1) / 2.0;
		y = RungeKutteMethod(alpha, max_error, false);
		if (y > B)
			alpha2 = alpha;
		else
			alpha1 = alpha;

		cout << setw(4) << fixed << setprecision(7)<< itr << " " << setw(12)<< alpha << " " << setw(12)<< y << " " << setprecision(10)<< max_error<< endl;

	} while (abs(y-B)>eps);
	RungeKutteMethod(alpha, max_error, true);
}

void main()
{
	shooting_metod();
}