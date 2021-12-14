#include<iostream>
#include<fstream>
#include<iomanip>
using namespace std;

const double eps1 = 0.00001; // для автоматического выбора шага
const double eps = 0.0001; // точность выполнения граничного условия

const int variant = 4;


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

void RungeKutteMethod()
{

}


void main()
{

}