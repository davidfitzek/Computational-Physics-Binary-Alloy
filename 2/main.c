#include <gsl/gsl_rng.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define E_cucu -0.436     // eV
#define E_znzn -0.113     // eV
#define E_cuzn -0.294     // eV
#define k_b 8.6173303e-05 // [eV/K]

gsl_rng *initialize_rng();
void equilibrate(gsl_rng *random_number, int nnList[][8], int *latticeA, int *latticeB, int N, int runs, double T);
double calc_energy(int nnList[][8], int *latticeA, int *latticeB, int N);
void init_bcc(int nnList[][8], int n, int N);
double calc_P(int *latticeA, int N);
double calc_mean(double *array, int length);
void write_vector(char *file_name, double *vec, int n, double dt, double t0);
void write_vector_int(char *file_name, int *vec, int n, double dt, double t0);
void metropolis(gsl_rng *random_number, int nnList[][8], int *latticeA, int *latticeB, int N, int runs, double T,
                int *accept_counter, double *P_iter, double *U_iter, double *r_iter);
double calc_r(int *latticeA, int *latticeB, int N, int nnList[][8]);
void calc_c(double *C, double *U, double delta_T, int T_len);
double calc_var(double *array, double mean, int N);
int calc_s_auto(double *array, int N);

int main()
{
    int n = 10;                                 // unit cells
    int N = n * n * n;                          // total number of cells
    int nnList[N][8];                           // nearest neighbour list
    int T;                                      // K
    int T_min = 0;                              // K
    int T_max = 1000;                           // K
    int delta_T = 50;                           // K
    int T_len = (int)(T_max - T_min) / delta_T; // K
    int init_iterations = 300000;               // number of iterations for equilibration
    int sim_iterations = 2000000;               // number of iterations for simulation
    int accept_counter;
    // allocate arrays
    int *latticeA = malloc(sizeof(double) * N);
    int *latticeB = malloc(sizeof(double) * N);
    double *P = malloc(sizeof(double) * T_len); // no unit
    double *P_iter = malloc(sizeof(double) * sim_iterations);
    double *U = malloc(sizeof(double) * T_len); // eV
    double *U_iter = malloc(sizeof(double) * sim_iterations);
    double *r = malloc(sizeof(double) * T_len); // no unit
    double *r_iter = malloc(sizeof(double) * sim_iterations);
    double *C = malloc(sizeof(double) * T_len - 1); // eV/K
    int *s_P = malloc(sizeof(double) * T_len);      // no unit

    gsl_rng *random_numbers = initialize_rng(); // setup random number generator
    init_bcc(nnList, n, N);                     // create neighbour matrix (N x 8)

    for (size_t iteration = 0; iteration < T_len; iteration++)
    {
        T = iteration * delta_T; // current temperature
        accept_counter = 0;
        // initialize 3D lattice to perfect order
        for (size_t i = 0; i < N; i++)
        {
            latticeA[i] = 0; // 0 = Cu
            latticeB[i] = 1; // 1 = Zn
        }

        // equlibrate system
        equilibrate(random_numbers, nnList, latticeA, latticeB, N, init_iterations, T);

        // run simulation
        metropolis(random_numbers, nnList, latticeA, latticeB, N, sim_iterations, T, &accept_counter, P_iter, U_iter, r_iter);

        P[iteration] = calc_mean(P_iter, sim_iterations);
        U[iteration] = calc_mean(U_iter, sim_iterations);
        r[iteration] = calc_mean(r_iter, sim_iterations);

        // compute correlation function
        s_P[iteration] = calc_s_auto(P_iter, sim_iterations);

        printf("P = %.4f, accepted = %.4f, r = %.4f, s_P = %.0d, T = %.4d\n", P[iteration], (double)accept_counter / sim_iterations,
               r[iteration], s_P[iteration], T);
    }

    calc_c(C, U, delta_T, T_len); // compute heat capacity finite difference approximation

    write_vector("data/P.dat", P, T_len, delta_T, T_min);
    write_vector("data/U.dat", U, T_len, delta_T, T_min);
    write_vector("data/r.dat", r, T_len, delta_T, T_min);
    write_vector("data/C.dat", C, T_len - 1, delta_T, T_min);

    write_vector_int("data/s_P.dat", s_P, T_len, delta_T, T_min);

    free(P);
    free(P_iter);
    free(U);
    free(U_iter);
    free(r);
    free(r_iter);
    free(latticeA);
    free(latticeB);
    free(C);

    return 0;
}

gsl_rng *initialize_rng()
{
    const gsl_rng_type *T;
    gsl_rng *q;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    q = gsl_rng_alloc(T);
    gsl_rng_set(q, time(NULL));
    return q;
}

void equilibrate(gsl_rng *random_number, int nnList[][8], int *latticeA, int *latticeB, int N, int runs, double T)
{
    double u, u_try, delta_u, q, r3;
    int r1, r2;
    int temp;
    int accept_counter = 0;

    // calc energy of system
    u = calc_energy(nnList, latticeA, latticeB, N);

    for (size_t i = 0; i < runs; i++)
    {
        // draw random numbers
        r1 = (int)floor(gsl_rng_uniform(random_number) * (N + 1));
        r2 = (int)floor(gsl_rng_uniform(random_number) * (N + 1));
        //int r = round(gsl_rng_uniform(random_number) * 999);
        r3 = gsl_rng_uniform(random_number);

        // change lattices
        temp = latticeA[r1];
        latticeA[r1] = latticeB[r2];
        latticeB[r2] = temp;

        // calc new energy of system
        u_try = calc_energy(nnList, latticeA, latticeB, N);

        delta_u = u - u_try;
        q = exp(delta_u / (k_b * T));
        if (q > r3)
        {
            accept_counter++;
            u = u_try;
        }
        else // switch back
        {
            temp = latticeA[r1];
            latticeA[r1] = latticeB[r2];
            latticeB[r2] = temp;
        }
    }
}

double calc_energy(int nnList[][8], int *latticeA, int *latticeB, int N)
{
    int index;
    double u = 0.0;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 8; j++)
        {
            index = nnList[i][j];
            if (latticeA[i] + latticeB[index] == 0)
                u += E_cucu;
            else if (latticeA[i] + latticeB[index] == 1)
                u += E_cuzn;
            else
                u += E_znzn;
        }
    }
    return u;
}

void init_bcc(int nnList[][8], int n, int N)
{
    for (size_t i = 0; i < N; i++)
    {
        nnList[i][0] = i;
        nnList[i][1] = (i + 1) % N;
        nnList[i][2] = (i + n) % N;
        nnList[i][3] = (i + n + 1) % N;
        nnList[i][4] = (i + n * n) % N;
        nnList[i][5] = (i + n * n + 1) % N;
        nnList[i][6] = (i + n * n + n) % N;
        nnList[i][7] = (i + n * n + n + 1) % N;
    }
}

double calc_P(int *latticeA, int N)
{
    double sum = 0;
    for (size_t i = 0; i < N; i++)
        sum += latticeA[i];
    //return 2 * (N - sum) / N - 1;
    return 1 - 2 * sum / N;
}

double calc_mean(double *array, int length)
{
    double sum = 0.0;
    for (size_t i = 0; i < length; i++)
    {
        sum += array[i];
    }
    return (sum / length);
}

void write_vector(char *file_name, double *vec, int n, double dt, double t0)
{
    FILE *fp;
    fp = fopen(file_name, "w");
    for (size_t i = 0; i < n; ++i)
    {
        fprintf(fp, "%.10f ", i * dt + t0);
        fprintf(fp, "%.10f ", vec[i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void write_vector_int(char *file_name, int *vec, int n, double dt, double t0)
{
    FILE *fp;
    fp = fopen(file_name, "w");
    for (size_t i = 0; i < n; ++i)
    {
        fprintf(fp, "%.10f ", i * dt + t0);
        fprintf(fp, "%d ", vec[i]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void metropolis(gsl_rng *random_number, int nnList[][8], int *latticeA, int *latticeB, int N, int runs, double T,
                int *accept_counter, double *P_iter, double *U_iter, double *r_iter)
{
    double u, u_try, delta_u, q, r3;
    int r1, r2;
    int temp;

    // calc energy of system
    u = calc_energy(nnList, latticeA, latticeB, N);

    for (size_t i = 0; i < runs; i++)
    {
        // draw random numbers
        r1 = (int)floor(gsl_rng_uniform(random_number) * (N + 1));
        r2 = (int)floor(gsl_rng_uniform(random_number) * (N + 1));
        r3 = gsl_rng_uniform(random_number);

        // change lattices
        temp = latticeA[r1];
        latticeA[r1] = latticeB[r2];
        latticeB[r2] = temp;

        // calc new energy of system
        u_try = calc_energy(nnList, latticeA, latticeB, N);

        delta_u = u - u_try;
        q = exp(delta_u / (k_b * T));
        if (q > r3)
        {
            *accept_counter = *accept_counter + 1;
            u = u_try;
        }
        else // switch back
        {
            temp = latticeA[r1];
            latticeA[r1] = latticeB[r2];
            latticeB[r2] = temp;
        }
        P_iter[i] = calc_P(latticeA, N);
        U_iter[i] = u;
        r_iter[i] = calc_r(latticeA, latticeB, N, nnList);
    }
}

double calc_r(int *latticeA, int *latticeB, int N, int nnList[][8])
{
    int index;
    double sum = 0.0;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < 8; j++)
        {
            index = nnList[i][j];
            if (latticeA[i] + latticeB[index] == 1)
                sum++;
        }
    }
    sum /= N;
    return 1.0 / 4.0 * (sum - 4);
}

void calc_c(double *C, double *U, double delta_T, int T_len)
{
    for (size_t i = 0; i < T_len - 1; i++)
    {
        C[i] = (U[i + 1] - U[i]) / delta_T;
    }
}

double calc_var(double *array, double mean, int N)
{
    double var = 0.0;
    for (size_t i = 0; i < N; i++)
    {
        var += (array[i] - mean) * (array[i] - mean);
    }
    var /= N;
    return var;
}

int calc_s_auto(double *array, int N)
{
    double mean = calc_mean(array, N);
    double var = calc_var(array, mean, N);
    double phi;
    int counter = 0;
    int corr_min = 0;
    int corr_max = 200000;
    int s = (corr_max + corr_min) / 2;

    while (abs(corr_max - corr_min) > 2 && counter < 100)
    {
        counter++;
        phi = 0.0;
        for (size_t i = 0; i < N - s; i++)
        {
            phi += (array[i] - mean) * (array[i + s] - mean);
        }
        phi /= (N - s);
        phi /= var;

        if (phi > 0.135)
        {
            corr_min = s;
        }
        else
        {
            corr_max = s;
        }
        s = (corr_max + corr_min) / 2;
    }
    return s;
}