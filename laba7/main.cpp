#define _USE_MATH_DEFINES
#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<cmath>

using namespace std;

const double eps1 = 0.00001; // ��� ��������������� ������ ����
const double eps = 0.0001; // �������� ���������� ���������� �������

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


double f(double x, double y, double z)
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
	return 0.1 * sin(M_PI * x) * variant + 0.1 * ksi * t * variant * pow(M_PI, 2) * sin(M_PI * x);
}

double FindH(double h, double alpha, double x, double y, double& z) /*� ����� ������� ��� ��� ��������*/
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

void FDM()
{
	for(int n=8;n<=32;n*=2)
	{
		double Del_T = 0;
		double h = 1. / n;
		double tau = tauExp(h);
		vector<double> cur_x(n+1);
		vector<double> next_x(n + 1);

		for (int i = 0; i < cur_x.size(); ++i)
			cur_x[i] = h * i;

		next_x = cur_x;

		for(int i=1;i*tauExp(h)<=1;++i)
		{
			//cur_x = ExplicitShem(n, i);
			double delta = 0;

			for (int j = 1; j < cur_x.size()-1; ++j)
			{
				next_x[j] = cur_x[j] + tau * ksi * (cur_x[j + 1] - 2 * cur_x[j] + cur_x[j - 1]) /(h*h) + tau * f2(tau*i, j * h);
			}

			cur_x= next_x;

			double max = 0;

			for (int j = 0; j <cur_x.size()-1; ++j)
			{
				if (abs(utx(tau*i, j * h) - cur_x[j]) > max)
					max = abs(utx(tau*i, j * h) - cur_x[j]);
				delta = max;
			}

			for(int i=0;i<cur_x.size();++i)
			{
				//cout << fixed << setprecision(5) << "|" << cur_x[i];
			}

			//cout << endl;

			if (delta > Del_T)
				Del_T = delta;

			//cout << tau*i << endl;
		}

		cout << Del_T << endl;
	}

	for (int n = 8; n <= 32; n *= 2)
	{
		double Del_T = 0;
		double h = 1. / n;
		double tau = tauImp(h);
		double d = tau * ksi / pow(h, 2);

		vector<double> cur_x(n + 1);
		vector<double> next_x(n + 1);

		for (int i = 0; i < cur_x.size(); ++i)
			cur_x[i] = h * i;

		next_x = cur_x;

		vector<vector<double>> a(n - 1, (vector<double>(n - 1)));

		for (int i = 0; i < a.size(); ++i)
		{
			if (i - 1 >= 0)
				a[i][i - 1] = -d;

			a[i][i] = 1 + 2 * d;

			if (i + 1 < a.size())
				a[i][i + 1] = -d;
		}

		for (int i = 1; i * tauImp(h) <= 1; ++i)
		{

			vector<double> b(n - 1);

			for (int j = 0; j < b.size(); ++j)
				b[j] = cur_x[j + 1] + tau * f2(tau * (i+1), h * (j+1) );

			b[b.size() - 1] += d;
			b = solveSLAU(a, b);

			for (int j = 0; j < b.size(); ++j)
				cur_x[j + 1] = b[j];
			

			double delta = 0;

			for (int j = 0; j <cur_x.size()-1; ++j)
			{
				double max = 0;
				if (abs(utx(tau*i, j * h) - cur_x[j]) > max)
					max = abs(utx(tau*i, j * h) - cur_x[j]);
				delta = max;
			}

			if (delta > Del_T)
				Del_T = delta;
		}

		cout << Del_T << endl;
	}
}

double YandZ(double h, double x, double y, double& z) /*���������� y � z*/
{
	double k11, k12, k13, k14; // y`=z
	double k21, k22, k23, k24; //z' = -Az + By - �sin(y) + F

	k11 = h * z;
	k21 = h * f(x, y, z);
	k12 = h * (z + k21 / 2.);
	k22 = h * f(x + h / 2., y + k11 / 2., z + k21 / 2.);
	k13 = h * (z + k22 / 2.);
	k23 = h * f(x + h / 2., y + k12 / 2., z + k22 / 2.);
	k14 = h * (z + k23);
	k24 = h * f(x + h, y + k13, z + k23);

	y = y + (k11 + 2 * k12 + 2 * k13 + k14) / 6.;
	z = z + (k21 + 2 * k22 + 2 * k23 + k24) / 6.;

	return y;
}


double RungeKutteMethod(double alpha, double& max_error, bool printT)
{

	/*y`(x) = z(x)
	z`(x)=f(x,y(x),z(x))
	y(0) = 1
	y(1) = 2
	z(0) = alpha*/


	double x = 0;
	double y = 1;
	double z =alpha;

	double y_prev, z_prev;

	if (printT)
		cout << "     x     |   y(x)   |   Ypr(x)  |   z(x)   |    Delta" << endl;

	while (x < 1) { // �� ���������� ��������, � ������ ������ ��������� �������� ��������� ����� ����������
		double h = 0.2;
		double y1, y2, y3;
		double z1;
		y_prev = y;
		z_prev = z;

		/*����� ������������ h*/
		do
		{
			z_prev = z;
			z1 = z;
			h = h / 2;

			y1 = YandZ(h, x, y_prev, z1); //y[j+1]=y[j] + SUM(p[i]k[i](h,y[j]))
			y2 = YandZ(h / 2., x, y_prev, z_prev); //y[j+1/2]=y[j] + SUM(p[i]k[i](h/2,y[j]))
			y3 = YandZ(h / 2.,  x, y2, z_prev); //y*[j+1]=y[j+1/2] + SUM(p[i]k[i](h/2,y[j+1/2]))
		} while (abs(y1 - y3) > eps1);

		double error = abs(Ynp(x) - y);
		if (error > max_error)
			max_error = error;

		if (printT)
			cout << fixed << setprecision(5) << setw(10) << x << setw(10) << y << setw(10) << Ynp(x) << setw(10) << z
			<< setw(14) << scientific << error << endl;


		y = YandZ(h, x, y, z);
		x += h;
	}

	return y_prev;
}


void shooting_metod()
{

	/*y`(x) = z(x)
	z`(x)=f(x,y(x),z(x))
	y(0) = 1
	z(0) = alpha*/

	double y;
	double alpha1 = 0, alpha2 = 3, alpha;
	double B = 2;
	int itr = 0;
	double max_error;
	cout << "  itr  |   z(0)     |    y(1)     |     Delta" << endl;
	/*��� ��������� ��� ������*/
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
		max_error = -1; // �� �������� �������� �� ������ ����� ����� ������ �������� ��� "��������"
	} while (abs(y-B)>eps);
	RungeKutteMethod(alpha, max_error, true);
}

void main()
{
	shooting_metod();
	FDM();
}