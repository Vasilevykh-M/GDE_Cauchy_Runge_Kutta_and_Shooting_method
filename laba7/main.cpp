#include<iostream>
#include<fstream>
#include<iomanip>
using namespace std;

const double eps1 = 0.00001; // ��� ��������������� ������ ����
const double eps = 0.0001; // �������� ���������� ���������� �������

const int variant = 4;


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


double f(double x, double y, double z) /*�� ��� ��� ���*/
{
	return -A(x) * z + B(x) * y - C(x) * sin(y) + F(x);
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


double RungeKutteMethod(double alpha, double& max_error, bool printTable)
{
	double x = 0;
	double y = 1;
	double z =alpha;

	double y_prev, z_prev;
	double k1, k2, k3, k4; // ��� ��������� y' = z
	double l1, l2, l3, l4; // ��� ��������� z' = -Az + By - �sin(y) + F

	if (printTable)
		cout << "     x     |   y(x)   |   Ypr(x)  |   z(x)   |    delta" << endl;

	while (x < 1) { // � ������� �� ���� ����� ����� ������� ����
		double h = 0.2;
		double y1, y2, y3;
		double z1;
		y_prev = y;
		z_prev = z;

		/*� ����� ��� ��� ��������*/
		do
		{
			y_prev = y;
			z_prev = z;
			z1 = z;

			h = h / 2;


			y1 = FindH(h, alpha, x, y_prev, z1);
			y2 = FindH(h / 2., alpha, x, y_prev, z_prev);
			y3 = FindH(h / 2., alpha, x, y2, z_prev);
		} while (abs(y1 - y3) > eps1);

		double error = abs(Ynp(x) - y);
		if (error > max_error)
			max_error = error;

		if (printTable)
			cout << fixed << setprecision(5) << setw(10) << x << setw(10) << y << setw(10) << Ynp(x) << setw(10) << z
			<< setw(14) << scientific << error << endl;


		y_prev = y;

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
	cout << "  itr  |   z(0)     |    y(1)     |     delta" << endl;
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

	} while (abs(y-B)>eps);
	RungeKutteMethod(alpha, max_error, true);
}

void main()
{
	shooting_metod();
}