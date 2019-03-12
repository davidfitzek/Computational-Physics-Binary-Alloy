#include <stdio.h>
#include <math.h>

#define E_cucu -0.436     // eV
#define E_znzn -0.113     // eV
#define E_cuzn -0.294     // eV
#define k_b 8.6173303e-05 // [eV/K]

double compute_f_prime(int N, double P, double delta_E, double T);
double compute_delta_E();
double compute_E0(int N);
double find_min_energy(double delta_E, double T, double N);
void print_array(double A[], int len);
void write_array(char *file_name, double **arr, int nrows, int ncols, double dt, double t0);
void write_vector(char *file_name, double *vec, int n, double dt, double t0);
double compute_u(double E0, double delta_E, double P, double N);
double compute_heat_capacity(double U1, double U2, double delta_T);

int main()
{
    int N = 100; // number of atoms
    int delta_T = 10;
    int T_min = 0;                              // K
    int T_max = 1000;                           // K
    double delta_E, E0;                         // eV
    int T_len = (int)(T_max - T_min) / delta_T; // K
    double P[T_len];
    double T[T_len];
    double U[T_len];
    double C[T_len - 1];

    delta_E = compute_delta_E();
    E0 = compute_E0(N);

    for (int i = 0; i < T_len; i++)
    {
        T[i] = i * delta_T; // 0k to 1000K
        P[i] = find_min_energy(delta_E, T[i], N);
        U[i] = compute_u(E0, delta_E, P[i], N);
    }

    for (int i = 1; i < T_len; i++)
    {
        C[i] = compute_heat_capacity(U[i], U[i + 1], delta_T);
    }

    write_vector("data/p.dat", P, T_len, delta_T, T_min);
    write_vector("data/u.dat", U, T_len, delta_T, T_min);
    write_vector("data/c.dat", C, T_len - 1, delta_T, T_min);

    printf("Phase transition occurs at T_c = %.4fK \n", 2 * delta_E / k_b);

    return 0;
}

double compute_f_prime(int N, double P, double delta_E, double T)
{
    double f_prime;
    f_prime = -4 * N * P * delta_E + N * T * k_b * (-log(1 - P) + log(1 + P));
    return f_prime;
}

double compute_delta_E()
{
    double delta_E;
    delta_E = E_znzn + E_cucu - 2 * E_cuzn;
    return delta_E;
}

double compute_E0(int N)
{
    double E0 = 2 * N * (E_cucu + E_znzn + 2 * E_cuzn);
    return E0;
}

double find_min_energy(double delta_E, double T, double N)
{
    double P_temp;
    double P_upper;
    double P_lower;
    double f_upper;
    double f_lower;
    double f_new;
    double P;
    int counter = 0;
    P_upper = 1.0;
    P_lower = 0.0;

    while (P_upper - P_lower > 0.001 && counter < 20)
    {
        counter++;
        f_upper = compute_f_prime(N, P_upper, delta_E, T);
        P_temp = (P_upper + P_lower) / 2;
        f_new = compute_f_prime(N, P_temp, delta_E, T);

        if (f_upper * f_new <= 0)
        {
            P_lower = P_temp;
        }
        else
        {
            P_upper = P_temp;
        }
    }
    return P_upper;
}

void print_array(double A[], int len)
{
    for (int row = 0; row < len; row++)
    {
        printf("%.2f     ", A[row]);
    }
    printf("\n\n");
}

void write_array(char *file_name, double **arr, int nrows, int ncols, double dt, double t0)
{

    FILE *fp;
    fp = fopen(file_name, "w");
    for (int i = 0; i < ncols; ++i)
    {
        fprintf(fp, "%.10f ", i * dt + t0);
        for (int j = 0; j < nrows; j++)
        {
            fprintf(fp, "%.10f ", arr[j][i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void write_vector(char *file_name, double *vec, int n, double dt, double t0)
{
    FILE *fp;
    fp = fopen(file_name, "w");
    for (int i = 0; i < n; ++i)
    {
        fprintf(fp, "%.10f ", i * dt + t0);
        fprintf(fp, "%.10f ", vec[i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

double compute_u(double E0, double delta_E, double P, double N)
{
    double U = E0 - 2 * N * P * delta_E;
    return U;
}

double compute_heat_capacity(double U1, double U2, double delta_T)
{
    double C = (U2 - U1) / delta_T;
    return C;
}