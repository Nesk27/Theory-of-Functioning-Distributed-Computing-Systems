#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS 1e-40

double fact(int n)
{
	return (n < 2) ? 1 : n * fact(n - 1);
}

//Формула расчета Pj при m = 1
double P_j(int j, int N, double u, double lambda)
{
	/*double first = pow((u / lambda), j);
	double second = 1.0 / fact(j);
	double third = 1.0 / exp(u / lambda);
	return (first * second * third);*/

	double sum = 0;
	for (int l = 0; l <= N; l++)
	{
		sum += (pow(u / lambda, l)) * (1.0 / fact(l));
	}

	sum = 1.0 / (sum);
	return (pow(u / lambda, j)) * (1.0 / fact(j)) * sum;
}

//Вспомогательный расчеты для вычислений Q
double pi_r(int r, int t, int i, double lambda)
{
	double first = pow(i * lambda * t, r) / fact(r);
	double second = exp(- 1 * i * lambda * t);
	return first * second;
}

int diff(int x)
{
	if (x >= 0) {
		return 1;
	} else
	{
		return 0;
	}
	return -1;
}

double u_l(int l, int t, int m, int N, int i, double u)
{
	double first = pow(u * t, l) / fact(l);
	double second = diff(N - i - m) * pow(m, l) * exp(-1 * m * u * t);
	double third = diff(m - N + i) * pow(N - i, l) * exp(-1 * (N - i) * u * t);
	

	return first * (second + third);
}

//Вычисление Q(t)
double Q_i(int i, int m, int N, int t, int n, double lambda, double u)
{
	int l = 0;
	double first = 0, second = 0, previous_step;
	
	do {
		previous_step = first;
		second = 0;
		for (int r = 0; r <= (i - n + l); r++)
			second += pi_r(r, t, i, lambda);
		first += (u_l(l, t, m, N, i, u) * second);

		l++;
	} while (fabs(first - previous_step) >= EPS);
	
	/*for (int l = 0; l < 100; l++)
	{
		second = 0;
		for (int r = 0; r <= (i - n + l); r++)
		{
			second += pi_r(r, t, i, lambda);
		}
		first += (u_l(l, t, m, N, i, u) * second);
	}*/
	//printf("l = %d\n", l);
	//printf("first in Q = %.20lf\n", first);

	
	
	return first;
}

double R(int n, int N, int m, int t, double u, double lambda)
{
	double sum = 0; double mul = 1; 
	double P, Q;
	for (int i = n; i <= N; i++)
	{
		P = P_j(i, N, u, lambda);
		Q = Q_i(i, m, N, t, n, lambda, u);
		sum += P * Q;
	}

	return sum;
}

double U(int t, int n, int m, int N, double u, double lambda) 
{
	double third = 0, second = 0, first = 0;
	for (int i = 0; i <= (n - 1); i++)
	{
		second = 0;
		for (int r = 0; r < 2; r++)
		{
			third = 0;
			for (int l = 0; l <= (n - i - 1 + r); l++)
			{
				third += u_l(l, t, m, N, i, u);
			}
			second += pi_r(r, t, i, lambda) * third;
		}
		first += P_j(i, N, u, lambda) * second;
	}

	return 1 - first;
}

double S(int N, int n, int m, double lambda, double u)
{
	/*double first = pow(lambda, (N - n + 1)) * pow((lambda + u), (-1 * (N - n + 1)));
	return 1 - first;*/
	/*double sum = 0, tmp = 0;
	for (int j = 0; j <= (n - 1); j++)
	{
		tmp = pow( (m * u) / lambda, j );
		sum += tmp * (1.0 / fact(j)) * exp(-1 * (m * u) / lambda);
	}
	return 1 - sum;*/
	double sum = 0;
	for (int j = n; j <= N; j++)
	{
		sum += P_j(j, N, u, lambda);
	}

	return sum;
}

int main()
{
	int N = 16, m = 16; //n = 10 /* n = [8..10] */, m = 1;
	double lambda = 0.024, u = 0.71; 
	//t = [0..24]

	//for(int j = 0; j <= N; j++)
	//printf("%d = %.20lf\n", j, P(j, u, lambda));

	/*for (int r = 0; r <= N; r++)
	{
		for (int i = 0; i < N; i++)
		{
			for (int t = 0; t <= 24; t++)
				printf("[r%d, i%d, t%d] = %.20lf\n", r, i, t, pi_r(r, t, i, lambda));
		}
		
	}*/

	/*for (int l = 0; l <= N; l++)
	{
		for (int i = n; i <= N; i++)
		{
			for (int t = 0; t <= 24; t++)
				printf("[l%d, i%d, t%d] = %.20lf\n", l, i, t, u_l(l, t, m, N, i, u));
		}

	}*/

	/*for (int i = n; i <= N; i++)
	{
		for (int t = 3; t <= 24; t++)
			printf("%.20lf\n", Q_i(i, m, N, t, n, lambda, u));
	}*/
	
	for (int n = 11; n <= 16; n++)
	{
		//printf("n = %d\n", n);
		printf("%.20lf\n", S(N, n, 16, lambda, u));

		//printf("[n%d, t%d] = %.20lf\n", n, t, R(n, N, m, t, u, lambda));
	}

	return EXIT_SUCCESS;
}