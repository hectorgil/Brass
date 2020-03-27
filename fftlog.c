#include <math.h>
#define Pi (4.*atan(1.))
#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>
#include "fftlog.h"
#include <stdio.h>

/* Computes the Gamma function using the Lanczos approximation */
double complex gammaFFT(double complex z) {
    /* Lanczos coefficients for g = 7 */
    double p[] = {
        0.99999999999980993227684700473478,
        676.520368121885098567009190444019,
       -1259.13921672240287047156078755283,
        771.3234287776530788486528258894,
       -176.61502916214059906584551354,
        12.507343278686904814458936853,
       -0.13857109526572011689554707,
        9.984369578019570859563e-6,
        1.50563273514931155834e-7
    };
    int n; 
    if(creal(z) < 0.5)
        return Pi / (csin(Pi*z)*gammaFFT(1. - z));
    z -= 1;
    double complex x = p[0];
    for(n = 1; n < 9; n++)
      x += p[n] / (z + (double)(n));
    double complex t = z + 7.5;
    return sqrt(2*Pi) * cpow(t, z+0.5) * cexp(-t) * x;
}

double complex polar (double r, double phi) {
  return (r*cos(phi) +I*(r*sin(phi)));
}

double complex lngamma(double complex z) {
    return clog(gammaFFT(z));
}


void lngamma_4(double x, double y, double* lnr, double* arg) {
    double complex w = lngamma(x+y*I);
    if(lnr) *lnr = creal(w);
    if(arg) *arg = cimag(w);
}

double goodkr(int N, double mu, double q, double L, double kr) {
    double xp = (mu+1+q)/2;
    double xm = (mu+1-q)/2;
    double y = Pi*N/(2*L);
    double lnr, argm, argp;
    lngamma_4(xp, y, &lnr, &argp);
    lngamma_4(xm, y, &lnr, &argm);
    double arg = log(2/kr) * N/L + (argp + argm)/Pi;
    double iarg = round(arg);
    if(arg != iarg)
        kr *= exp((arg - iarg)*L/N);
    return kr;
}

void compute_u_coefficients(int N, double mu, double q, double L, double kcrc, double complex u[]) {
    double y = Pi/L;
    double k0r0 = kcrc * exp(-L);
    double t = -2*y*log(k0r0/2);
    int m;
    if(q == 0) {
        double x = (mu+1)/2;
        double lnr, phi;
        for( m = 0; m <= N/2; m++) {
            lngamma_4(x, m*y, &lnr, &phi);
            u[m] = polar(1.0,m*t + 2*phi);
        }
    }
    else {
        double xp = (mu+1+q)/2;
        double xm = (mu+1-q)/2;
        double lnrp, phip, lnrm, phim;
        for( m = 0; m <= N/2; m++) {
            lngamma_4(xp, m*y, &lnrp, &phip);
            lngamma_4(xm, m*y, &lnrm, &phim);
            u[m] = polar(exp(q*log(2) + lnrp - lnrm), m*t + phip - phim);
        }
    }

    for( m = N/2+1; m < N; m++)
        u[m] = conj(u[N-m]);
    if((N % 2) == 0)
      u[N/2] = (creal(u[N/2]) + I*0.0);
}

void fht_threadsafe(int N, double r[], double complex a[], double k[], double complex b[], double mu, double q, double kcrc, int noring, double complex* u,fftw_plan p_forward, fftw_plan p_reverse)
{
    int m,n;
    double L = log(r[N-1]/r[0]) * N/(N-1.);
    double complex* ulocal = NULL;
    if(u == NULL) {
        if(noring)
            kcrc = goodkr(N, mu, q, L, kcrc);
        ulocal = malloc (sizeof(complex double)*N);
        compute_u_coefficients(N, mu, q, L, kcrc, ulocal);
        u = ulocal;
    }
      fftw_execute_dft(p_forward,a,b);
    for( m = 0; m < N; m++)
      b[m] *= u[m] / (double)(N);       // divide by N since FFTW doesn't normalize the inverse FFT
      fftw_execute_dft(p_reverse,b,b);

     

    double complex tmp;
    for( n = 0; n < N/2; n++) {
        tmp = b[n];
        b[n] = b[N-n-1];
        b[N-n-1] = tmp;
    }

    double k0r0 = kcrc * exp(-L);
    k[0] = k0r0/r[0];
    for( n = 1; n < N; n++)
        k[n] = k[0] * exp(n*L/N);

    free(ulocal);
}


void fht(int N, double r[], double complex a[], double k[], double complex b[], double mu, double q, double kcrc, int noring, double complex* u)
{
    int m,n;
    double L = log(r[N-1]/r[0]) * N/(N-1.);
    double complex* ulocal = NULL;
    if(u == NULL) {
        if(noring)
            kcrc = goodkr(N, mu, q, L, kcrc);
        ulocal = malloc (sizeof(complex double)*N); 
        compute_u_coefficients(N, mu, q, L, kcrc, ulocal);
        u = ulocal;
    }

    /* Compute the convolution b = a*u using FFTs */
    fftw_plan forward_plan = fftw_plan_dft_1d(N, (fftw_complex*) a, (fftw_complex*) b,  -1, FFTW_ESTIMATE);
    fftw_plan reverse_plan = fftw_plan_dft_1d(N, (fftw_complex*) b, (fftw_complex*) b, +1, FFTW_ESTIMATE);
    fftw_execute(forward_plan);
       
    for( m = 0; m < N; m++)
      b[m] *= u[m] / (double)(N);       // divide by N since FFTW doesn't normalize the inverse FFT
    fftw_execute(reverse_plan);

    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(reverse_plan);

    /* Reverse b array */
    double complex tmp;
    for( n = 0; n < N/2; n++) {
        tmp = b[n];
        b[n] = b[N-n-1];
        b[N-n-1] = tmp;
    }

    /* Compute k's corresponding to input r's */
    double k0r0 = kcrc * exp(-L);
    k[0] = k0r0/r[0];
    for( n = 1; n < N; n++)
        k[n] = k[0] * exp(n*L/N);

    free(ulocal);
    fftw_cleanup();
}

/* Compute the function  \xi_l^m(r) = \int_0^\infty \frac{dk}{2\pi^2} k^m j_l(kr) P(k)
    Note that the usual 2-point correlation function xi(r) is just xi_0^2(r)
    in this notation.  The input k-values must be logarithmically spaced.  The
    resulting xi_l^m(r) will be evaluated at the dual r-values
    r[0] = 1/k[N-1], ..., r[N-1] = 1/k[0]. 
*/
void fftlog_ComputeXiLM(int l, int m, int N, double k[], double pk[], double r[], double xi[]) {
  double complex* a = malloc(sizeof(complex double)*N);//move to calloc?
  double complex* b = malloc(sizeof(complex double)*N);//move to calloc?
  int i;
    for( i = 0; i < N; i++)
        a[i] = pow(k[i], m - 0.5) * pk[i];//assigned real part only
    fht(N, k, a, r, b, l + 0.5, 0, 1, 1, NULL);//b modified, a modified?
    for( i = 0; i < N; i++)
        xi[i] = creal(pow(2*Pi*r[i], -1.5) * b[i]);

    free(a);
    free(b);
}

void fftlog_ComputeXiLM_threadsafe(int l, int m, int N, double k[], double pk[], double r[], double xi[],fftw_plan p_forward, fftw_plan p_reverse) {
  double complex* a = malloc(sizeof(complex double)*N);//move to calloc?
  double complex* b = malloc(sizeof(complex double)*N);//move to calloc?
  int i;
    for( i = 0; i < N; i++)
        a[i] = pow(k[i], m - 0.5) * pk[i];//assigned real part only
    fht_threadsafe(N, k, a, r, b, l + 0.5, 0, 1, 1, NULL,p_forward,p_reverse);//b modified, a modified?
    for( i = 0; i < N; i++)
        xi[i] = creal(pow(2*Pi*r[i], -1.5) * b[i]);

    free(a);
    free(b);
}
