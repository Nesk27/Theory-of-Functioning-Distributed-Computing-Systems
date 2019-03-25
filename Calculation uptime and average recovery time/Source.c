#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>

int calc_ul(int N, int m, int u, int l)
{
	if (((N - m) <= l) && (l <= N))
	{
		return (N - l) * u;
	} else {
		return m * u;

	};


	/*if (l < (N - m))
		return m * u;

	return (N - l) * u;*/

	
}

double calc_T(double lambda, int N, int n, int m, int u)
{
	double T, iter1 = 1, iter2 = 0, iter3 = 1;
	//—читаем первое произведение в формуле
	for (int l = 1; l <= (n - 1); l++)
		iter1 *= (l * lambda) / calc_ul(N, m, u, l);

	iter1 /= calc_ul(N, m, u, 0);

	//—читаем второе слагаемое в формуле (¬нешн€€ сумма)
	for (int j = 1; j <= (n - 1); j++)
	{
		iter3 = 1;
		for (int l = j; l <= (n - 1); l++)
		{
			iter3 *= (l * lambda) / calc_ul(N, m, u, l);
		}
		iter2 += (1 / (j * lambda)) * iter3;
	}

	T = iter1 + iter2;
	return T;
}

double calc_Theta(double lambda, int N, int n, int m, int u)
{
	double iter1 = 0, iter2 = 1, iter3 = 0;
	for (int j = (n + 1); j <= N; j++)
	{
		
		
		iter1 = 1.0 / (j * lambda);
		iter2 = 1;
		for (int l = n; l <= (j - 1); l++)
		{
			iter2 *= calc_ul(N, m, u, l) / (l * lambda);
		}
		iter3 += iter1 * iter2;
	}
	return iter3 + 1.0 / (n * lambda);

}


int main()
{
	/*double N = 8192, m = 4, u = 1;
	double lambda = 0.00001;*/
	double T;
		//T = calc_T(0.00001, 8192, 8192, 1, 1);
		//printf("%lf\n", T);
	for (int m = 1; m <= 4; m++) {
		printf("m = %d\n", m);
		for (int n = 65527; n <= 65536; n++)
		{
			double Theta = calc_Theta(0.00001, 65536, n, m, 1);
			printf("%f\n", Theta);
		}
	}

	
	return 0;

}