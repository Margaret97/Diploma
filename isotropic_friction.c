#include "pch.h"
#include <iostream>
#include <stdlib.h> 
#define _USE_MATH_DEFINES 
#include <math.h> 
#include<float.h>
#define delta_1 1e-11
#define delta_2 1e-11
#define accuracy 1e-13
#define eps 1e-11

double f_1(double y); // the right side of the first equation
double f_2(double x, double y, double p, double epsilon); // the right side of the second equation
void Runge_Kutta_8(double* x, double* y, double step, double p, double epsilon); // the Runge-Kutta method of the 8th-order
double local_error(double actual_x, double actual_y, double step, double p, double epsilon);// the local error for the Runge-Kutta method
double max(double a, double b); // the maximum function
double min(double a, double b); // the minimum function
int find_cycle_period(double p, double epsilon, double* T); // search of the period

double f_1(double y)
{

	return y;

}
double f_2(double x, double y, double p, double epsilon) // p corresponds to k in diploma
{

	return -p * (y - 1) - epsilon * p * sin(x);

}
void Runge_Kutta_8(double* x, double* y, double step, double p, double epsilon)
{
	double k[13], q[13];
	k[0] = step * f_1(*y);
	q[0] = step * f_2(*x, *y, p, epsilon);
	k[1] = step * f_1(*y + (1 / 18.) * q[0]);
	q[1] = step * f_2(*x + (1 / 18.) * k[0], *y + (1 / 18.) * q[0], p, epsilon);
	k[2] = step * f_1(*y + (1 / 48.) * q[0] + (1 / 16.) * q[1]);
	q[2] = step * f_2(*x + (1 / 48.) * k[0] + (1 / 16.) * k[1], *y + (1 / 48.) * q[0] + (1 / 16.) * q[1], p, epsilon);
	k[3] = step * f_1(*y + (1 / 32.) * q[0] + (3 / 32.) * q[2]);
	q[3] = step * f_2(*x + (1 / 32.) * k[0] + (3 / 32.) * k[2], *y + (1 / 32.) * q[0] + (3 / 32.) * q[2], p, epsilon);
	k[4] = step * f_1(*y + (5 / 16.) * q[0] - (75 / 64.) * q[2] + (75 / 64.) * q[3]);
	q[4] = step * f_2(*x + (5 / 16.) * k[0] - (75 / 64.) * k[2] + (75 / 64.) * k[3], *y + (5 / 16.) * q[0] - (75 / 64.) * q[2] + (75 / 64.) * q[3], p, epsilon);
	k[5] = step * f_1(*y + (3 / 80.) * q[0] + (3 / 16.) * q[3] + (3 / 20.) * q[4]);
	q[5] = step * f_2(*x + (3 / 80.) * k[0] + (3 / 16.) * k[3] + (3 / 20.) * k[4], *y + (3 / 80.) * q[0] + (3 / 16.) * q[3] + (3 / 20.) * q[4], p, epsilon);
	k[6] = step * f_1(*y + (29443841 / 614563906.) * q[0] + (77736538 / 692538347.) * q[3] - (28693883 / 1125000000.) * q[4] + (23124283 / 1800000000.) * q[5]);
	q[6] = step * f_2(*x + (29443841 / 614563906.) * k[0] + (77736538 / 692538347.) * k[3] - (28693883 / 1125000000.) * k[4] + (23124283 / 1800000000.) * k[5], *y + (29443841 / 614563906.) * q[0] + (77736538 / 692538347.) * q[3] - (28693883 / 1125000000.) * q[4] + (23124283 / 1800000000.) * q[5], p, epsilon);
	k[7] = step * f_1(*y + (16016141 / 946692911.) * q[0] + (61564180 / 158732637.) * q[3] + (22789713 / 633445777.) * q[4] + (545815736 / 2771057229.) * q[5] - (180193667 / 1043307555.) * q[6]);
	q[7] = step * f_2(*x + (16016141 / 946692911.) * k[0] + (61564180 / 158732637.) * k[3] + (22789713 / 633445777.) * k[4] + (545815736 / 2771057229.) * k[5] - (180193667 / 1043307555.) * k[6], *y + (16016141 / 946692911.) * q[0] + (61564180 / 158732637.) * q[3] + (22789713 / 633445777.) * q[4] + (545815736 / 2771057229.) * q[5] - (180193667 / 1043307555.) * q[6], p, epsilon);
	k[8] = step * f_1(*y + (39632708 / 573591083.) * q[0] - (433636366 / 683701615.) * q[3] - (421739975 / 2616292301.) * q[4] + (100302831 / 723423059.) * q[5] + (790204164 / 839813087.) * q[6] + (800635310 / 3783071287.) * q[7]);
	q[8] = step * f_2(*x + (39632708 / 573591083.) * k[0] - (433636366 / 683701615.) * k[3] - (421739975 / 2616292301.) * k[4] + (100302831 / 723423059.) * k[5] + (790204164 / 839813087.) * k[6] + (800635310 / 3783071287.) * k[7], *y + (39632708 / 573591083.) * q[0] - (433636366 / 683701615.) * q[3] - (421739975 / 2616292301.) * q[4] + (100302831 / 723423059.) * q[5] + (790204164 / 839813087.) * q[6] + (800635310 / 3783071287.) * q[7], p, epsilon);
	k[9] = step * f_1(*y + (246121993 / 1340847787.) * q[0] - (37695042795 / 15268766246.) * q[3] - (309121744 / 1061227803.) * q[4] - (12992083 / 490766935.) * q[5] + (6005943493 / 2108947869.) * q[6] + (393006217 / 1396673457.) * q[7] + (123872331 / 1001029789.) * q[8]);
	q[9] = step * f_2(*x + (246121993 / 1340847787.) * k[0] - (37695042795 / 15268766246.) * k[3] - (309121744 / 1061227803.) * k[4] - (12992083 / 490766935.) * k[5] + (6005943493 / 2108947869.) * k[6] + (393006217 / 1396673457.) * k[7] + (123872331 / 1001029789.) * k[8], *y + (246121993 / 1340847787.) * q[0] - (37695042795 / 15268766246.) * q[3] - (309121744 / 1061227803.) * q[4] - (12992083 / 490766935.) * q[5] + (6005943493 / 2108947869.) * q[6] + (393006217 / 1396673457.) * q[7] + (123872331 / 1001029789.) * q[8], p, epsilon);
	k[10] = step * f_1(*y - (1028468189 / 846180014.) * q[0] + (8478235783 / 508512852.) * q[3] + (1311729495 / 1432422823.) * q[4] - (10304129995 / 1701304382.) * q[5] - (48777925059 / 3047939560.) * q[6] + (15336726248 / 1032824649.) * q[7] - (45442868181 / 3398467696.) * q[8] + (3065993473 / 597172653.) * q[9]);
	q[10] = step * f_2(*x - (1028468189 / 846180014.) * k[0] + (8478235783 / 508512852.) * k[3] + (1311729495 / 1432422823.) * k[4] - (10304129995 / 1701304382.) * k[5] - (48777925059 / 3047939560.) * k[6] + (15336726248 / 1032824649.) * k[7] - (45442868181 / 3398467696.) * k[8] + (3065993473 / 597172653.) * k[9], *y - (1028468189 / 846180014.) * q[0] + (8478235783 / 508512852.) * q[3] + (1311729495 / 1432422823.) * q[4] - (10304129995 / 1701304382.) * q[5] - (48777925059 / 3047939560.) * q[6] + (15336726248 / 1032824649.) * q[7] - (45442868181 / 3398467696.) * q[8] + (3065993473 / 597172653.) * q[9], p, epsilon);
	k[11] = step * f_1(*y + (185892177 / 718116043.) * q[0] - (3185094517 / 667107341.) * q[3] - (477755414 / 1098053517.) * q[4] - (703635378 / 230739211.) * q[5] + (5731566787 / 1027545527.) * q[6] + (5232866602 / 850066563.) * q[7] - (4093664535 / 808688257.) * q[8] + (3962137247 / 1805957418.) * q[9] + (65686358 / 487910083.) * q[10]);
	q[11] = step * f_2(*x + (185892177 / 718116043.) * k[0] - (3185094517 / 667107341.) * k[3] - (477755414 / 1098053517.) * k[4] - (703635378 / 230739211.) * k[5] + (5731566787 / 1027545527.) * k[6] + (5232866602 / 850066563.) * k[7] - (4093664535 / 808688257.) * k[8] + (3962137247 / 1805957418.) * k[9] + (65686358 / 487910083.) * k[10], *y + (185892177 / 718116043.) * q[0] - (3185094517 / 667107341.) * q[3] - (477755414 / 1098053517.) * q[4] - (703635378 / 230739211.) * q[5] + (5731566787 / 1027545527.) * q[6] + (5232866602 / 850066563.) * q[7] - (4093664535 / 808688257.) * q[8] + (3962137247 / 1805957418.) * q[9] + (65686358 / 487910083.) * q[10], p, epsilon);
	k[12] = step * f_1(*y + (403863854 / 491063109.) * q[0] - (5068492393 / 434740067.) * q[3] - (411421997 / 543043805.) * q[4] + (652783627 / 914296604.) * q[5] + (11173962825 / 925320556.) * q[6] - (13158990841 / 6184727034.) * q[7] + (3936647629 / 1978049680.) * q[8] - (160528059 / 685178525.) * q[9] + (248638103 / 1413531060.) * q[10]);
	q[12] = step * f_2(*x + (403863854 / 491063109.) * k[0] - (5068492393 / 434740067.) * k[3] - (411421997 / 543043805.) * k[4] + (652783627 / 914296604.) * k[5] + (11173962825 / 925320556.) * k[6] - (13158990841 / 6184727034.) * k[7] + (3936647629 / 1978049680.) * k[8] - (160528059 / 685178525.) * k[9] + (248638103 / 1413531060.) * k[10], *y + (403863854 / 491063109.) * q[0] - (5068492393 / 434740067.) * q[3] - (411421997 / 543043805.) * q[4] + (652783627 / 914296604.) * q[5] + (11173962825 / 925320556.) * q[6] - (13158990841 / 6184727034.) * q[7] + (3936647629 / 1978049680.) * q[8] - (160528059 / 685178525.) * q[9] + (248638103 / 1413531060.) * q[10], p, epsilon);
	*x = *x + (14005451 / 335480064.) * k[0] + (-59238493 / 1068277825.) * k[5] + (181606767 / 758867731.) * k[6] + (561292985 / 797845732.) * k[7] - (1041891430 / 1371343529.) * k[8] + (760417239 / 1151165299.) * k[9] + (118820643 / 751138087.) * k[10] - (528747749 / 2220607170.) * k[11] + (1 / 4.) * k[12];
	*y = *y + (14005451 / 335480064.) * q[0] + (-59238493 / 1068277825.) * q[5] + (181606767 / 758867731.) * q[6] + (561292985 / 797845732.) * q[7] - (1041891430 / 1371343529.) * q[8] + (760417239 / 1151165299.) * q[9] + (118820643 / 751138087.) * q[10] - (528747749 / 2220607170.) * q[11] + (1 / 4.) * q[12];
}

double local_error(double actual_x, double actual_y, double step, double p, double epsilon)
{
	double value_x[2], value_y[2];
	double error_x, error_y, error;
	value_x[0] = actual_x;
	value_y[0] = actual_y;
	Runge_Kutta_8(&value_x[0], &value_y[0], step * (1 / 2.), p, epsilon);
	Runge_Kutta_8(&value_x[0], &value_y[0], step * (1 / 2.), p, epsilon);
	value_x[1] = actual_x;
	value_y[1] = actual_y;
	Runge_Kutta_8(&value_x[1], &value_y[1], step, p, epsilon);
	error_x = (fabs(value_x[1] - value_x[0])) * (1 / 255.);
	error_y = (fabs(value_y[1] - value_y[0])) * (1 / 255.);
	error = max(error_x, error_y);
	return error;
}

double max(double a, double b)
{
	if (a > b)
		return a;
	else return b;
}

double min(double a, double b)
{
	if (a < b)
		return a;
	else return b;
}
int find_cycle_period(double p, double epsilon, double* T)
{
	FILE* t_x_y;
	double x, y, t, prev_t, prev_x, prev_y, x_0, y_0, t_0, current_t, x_add, y_add;
	double h; // the integration step
	double err;
	int i, j;
	double prev_T;
	int k;
	k = 0;
	t = 0; i = 1; j = 0;
	h = 0.01;
	x = 0; y = 1;
	prev_t = 0;
	*T = 2 * M_PI;
	fopen_s(&t_x_y, "C://Users//Margarita//data//task_2.txt", "w");

	do // We find the periodic solution (x and y:  y(x)=y(x+2\pi))
	{
		x_0 = x; y_0 = y;
		prev_x = x;
		while ((prev_x - x_0 - 2 * M_PI) * (x - x_0 - 2 * M_PI) > 0) // we integrate equations until the expression x-(x_0+2\pi) changes sign
		{
			prev_t = t;
			prev_x = x;
			prev_y = y;
			err = local_error(x, y, h, p, epsilon);
			if (err <= eps)
			{
				Runge_Kutta_8(&x, &y, 0.5 * h, p, epsilon);
				Runge_Kutta_8(&x, &y, 0.5 * h, p, epsilon);
				t = t + h;

				if (err < 1e-30)
				{
					h = 2 * h;
				}

				else
				{
					h = h * min(2, max((1 / 3.), 0.8 * pow((eps / err), (1 / 9.))));
				}
			}
			else
			{
				while (err > eps)
				{
					h = h * min(2, max((1 / 3.), 0.8 * pow((eps / err), (1 / 9.))));
					err = local_error(x, y, h, p, epsilon);
					i++;
				}
				Runge_Kutta_8(&x, &y, 0.5 * h, p, epsilon);
				Runge_Kutta_8(&x, &y, 0.5 * h, p, epsilon);
				t = t + h;
				if (err < 1e-30)
				{
					h = 2 * h;
				}

				else
				{
					h = h * min(2, max((1 / 3.), 0.8 * pow((eps / err), (1 / 9.))));
				}
			}


		}

		k = 0;
		// the golden-section search (it is less efficient than the secant method)
		/*while (fabs(prev_x - x_0 - 2 * M_PI) > accuracy && fabs(x - x_0 - 2 * M_PI) > accuracy)
		{
			x_add = prev_x; y_add = prev_y;
			current_t = prev_t + (t - prev_t) / 1.618;
			if (fabs(current_t - prev_t) < 1e-35)
			{
				printf("(cycle)Very small accuracy\n");
				printf("Amount of operation is %d\n", k);
				return -1;

			}
			Runge_Kutta_8(&prev_x, &prev_y, (current_t - prev_t), p, epsilon);
			if ((prev_x - x_0 - 2 * M_PI)*(x - x_0 - 2 * M_PI) > 0)
			{
				t = current_t;
				x = prev_x;
				y = prev_y;
				prev_x = x_add;
				prev_y = y_add;
			}
			else
			{
				prev_t = current_t;
			}
			k++;
		}*/
		current_t = prev_t;

		while (fabs(prev_x - x_0 - 2 * M_PI) > accuracy && fabs(x - x_0 - 2 * M_PI) > accuracy) // the secant method
		{
			x_add = prev_x; y_add = prev_y;
			current_t = prev_t - (prev_x - 2 * M_PI - x_0) * (t - prev_t) / (x - prev_x);
			if (fabs(current_t - prev_t) < 1e-35)
			{
				printf("(cycle)Very small accuracy\n");
				printf("Amount of operation is %d\n", k);
				return -1;

			}
			Runge_Kutta_8(&prev_x, &prev_y, (current_t - prev_t), p, epsilon);
			if ((prev_x - x_0 - 2 * M_PI) * (x - x_0 - 2 * M_PI) > 0)
			{
				t = current_t;
				x = prev_x;
				y = prev_y;
				prev_x = x_add;
				prev_y = y_add;

			}
			else
			{
				prev_t = current_t;
			}
			k++;


		}
		if (fabs(prev_x - x_0 - 2 * M_PI) <= accuracy)
		{
			t = current_t; x = prev_x; y = prev_y;
		}

		j++;
		printf("Step %d amount of operation is %d\n", j, k);
		printf("x_0+2pi*%d=%.16e\n", j, x);
		fprintf_s(t_x_y, "(cycle)Step %d\t%d\n", j, k);
		fprintf_s(t_x_y, "&$%.16e$\t&$%.16e$\n", t, x);

	} while (fabs(y - y_0) > delta_1);

	//printf("We will find period\n");


	j = 0;
	do // We find the period of the periodic solution (x and y:  y(x) = y(x + 2\pi))
	{
		x_0 = x;
		t_0 = t;
		prev_x = x;
		prev_T = *T;
		while ((prev_x - x_0 - 2 * M_PI) * (x - x_0 - 2 * M_PI) > 0) // we integrate equations until the expression x-(x_0+2\pi) changes sign
		{
			prev_t = t;
			prev_x = x;
			prev_y = y;
			err = local_error(x, y, h, p, epsilon);
			if (err <= eps)
			{
				Runge_Kutta_8(&x, &y, 0.5 * h, p, epsilon);
				Runge_Kutta_8(&x, &y, 0.5 * h, p, epsilon);
				t = t + h;
				if (err < 1e-30)
				{
					h = 2 * h;
				}

				else
				{
					h = h * min(2, max((1 / 3.), 0.8 * pow((eps / err), (1 / 9.))));
				}
			}
			else
			{
				i = 1;
				while (err > eps)
				{
					h = h * min(2, max((1 / 3.), 0.8 * pow((eps / err), (1 / 9.))));
					err = local_error(x, y, h, p, epsilon);
					i++;
				}
				Runge_Kutta_8(&x, &y, 0.5 * h, p, epsilon);
				Runge_Kutta_8(&x, &y, 0.5 * h, p, epsilon);
				t = t + h;
				if (err < 1e-30)
				{
					h = 2 * h;
				}

				else
				{
					h = h * min(2, max((1 / 3.), 0.8 * pow((eps / err), (1 / 9.))));
				}
			}


		}

		k = 0;
		// the golden-section search (it is less efficient than the secant method)
		/*while (fabs(prev_x - x_0 - 2 * M_PI) > accuracy && fabs(x - x_0 - 2 * M_PI) > accuracy)
		{
			x_add = prev_x; y_add = prev_y;
			current_t = prev_t + (t - prev_t) / 1.618;
			if (fabs(current_t - prev_t) < 1e-35)
			{
				printf("(period)Very small accuracy\n");
				printf("Amount of operation is %d\n", k);
				return -1;

			}
			Runge_Kutta_8(&prev_x, &prev_y, (current_t - prev_t), p, epsilon);
			if ((prev_x - x_0 - 2 * M_PI)*(x - x_0 - 2 * M_PI) > 0)
			{
				t = current_t;
				x = prev_x;
				y = prev_y;
				prev_x = x_add;
				prev_y = y_add;
			}
			else
			{
				prev_t = current_t;
			}
			k++;
		}*/
		current_t = prev_t;
		while (fabs(prev_x - x_0 - 2 * M_PI) > accuracy && fabs(x - x_0 - 2 * M_PI) > accuracy) // the secant method
		{
			x_add = prev_x; y_add = prev_y;
			current_t = prev_t - (prev_x - 2 * M_PI - x_0) * (t - prev_t) / (x - prev_x);
			if (fabs(current_t - prev_t) < 1e-30)
			{
				printf("(period)Very small accuracy\n");
				printf("Amount of operation is %d\n", k);
				return -1;

			}
			Runge_Kutta_8(&prev_x, &prev_y, (current_t - prev_t), p, epsilon);

			if ((prev_x - x_0 - 2 * M_PI) * (x - x_0 - 2 * M_PI) > 0)
			{
				t = current_t;
				x = prev_x;
				y = prev_y;
				prev_x = x_add;
				prev_y = y_add;
			}
			else
			{
				prev_t = current_t;
			}
			k++;


		}
		if (fabs(prev_x - x_0 - 2 * M_PI) <= accuracy)
		{
			t = current_t; x = prev_x; y = prev_y;
		}
		*T = t - t_0;
		j++;

		printf("Period: Step %d amount of operation is %d\n", j, k);
		fprintf_s(t_x_y, "(period)Step %d\t%d \n", j, k);
		fprintf_s(t_x_y, "&$%.16e$\t&$%.16e$\n", t, x);

	} while (fabs(*T - prev_T) > delta_2);
	fprintf_s(t_x_y, "&$%.16e$\n", *T);
	fclose(t_x_y);
	return 1;
}




int main(void)
{
	FILE* epsilon_period;
	int check;
	double p;
	double epsilon;
	double period;
	epsilon = 0.1;
	printf("Enter p:\n");
	scanf_s("%lf", &p);
	fopen_s(&epsilon_period, "C://Users//Margarita//data//epsilon_period.txt", "w");
	if (epsilon_period == NULL)
	{
		printf("Error opening file epsilon_period\n");
		return -1;
	}

	check = find_cycle_period(p, epsilon, &period);
	if (check == 1)
	{
		fprintf_s(epsilon_period, "%lf\t %lf\n", epsilon, period);
	}
	else
	{
		return -1;
	}

	fclose(epsilon_period);
	return 0;
}
