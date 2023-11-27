// gcc -Wall -I/usr/local/include/gsl -c  prog_06_01.cpp
// gcc -L/usr/local/lib prog_06_01.o -lgsl -lgslcblas -lm
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


// model parameters
double gamma_p = 0.5;
double egam = 1.0-1.0/(double)gamma_p;
double beta = 0.9;
double alpha = 0.3;
double delta = 0.0;
double by = 0.0;
double kappa = 0.0;
double n_p = 0.2;
signed short int tax = 1;
double g[] = {0.12, 0.12, 0.0};

// model variables

double w, r, wn, Rn, p, tauw, taur, tauc, taup, pen;
double KK, LL, YY, AA, CC, BB, GG, II;
double a[3];
double c[3];

int eqns_f(const gsl_vector * x,void *params, gsl_vector * f)
{
    double PVI, PSI;

    const double x0 = gsl_vector_get (x, 0);
    const double x1 = gsl_vector_get (x, 1);

    /*
    const double y0 = a * (1 - x0);
    const double y1 = b * (x1 - x0 * x0);
    */

    // copy values to respective variables
    KK = x0;

    if(tax == 1)
    {
        tauc = x1;
    }
    if(tax == 2)
    {
        tauw = x1;
        taur = x1;
    }
    if(tax == 3)
    {
        tauw = x1;
    }
    if(tax > 3)
    {
        taur = x1;
    }

    // factor prices and pension payments
    r = alpha*pow(KK/LL,(alpha-1.0))-delta;
    w = (1.0-alpha)*pow(KK/LL,alpha);
    wn = w*(1.0-tauw-taup);
    Rn = 1.0+r*(1.0-taur);
    p = 1.0+tauc;
    pen = kappa*w;

    // individual decisions
    PVI = wn + wn/Rn + pen/pow(Rn,2);
    PSI = p*(1.0 + pow(beta*Rn,gamma_p)/Rn + pow(beta*Rn,2.0*gamma_p)/pow(Rn,2));
    c[0] = PVI/PSI;
    c[1] = pow((beta*Rn),gamma_p)*c[0];
    c[2] = pow((beta*Rn),gamma_p)*c[1];
    a[0] = 0.0;
    a[1] = wn - p*c[0];
    a[2] = wn + Rn*a[1] - p*c[1];

    // quantities
    YY = pow(KK,alpha) * pow(LL,1-0-alpha);
    CC = c[0] + c[1]/(1.0+n_p) + c[2]/pow(1.0+n_p,2);
    GG = g[0] + g[1]/(1.0+n_p) + g[2]/pow(1.0+n_p,2);
    AA = a[1]/(1.0+n_p) + a[2]/pow(1.0+n_p,2);
    II = (n_p+delta)*KK;
    BB = by*YY;

    // get equations defining general equilibrium
    const double y0 = KK + BB - AA;
    const double y1 = tauc*CC + tauw*w*LL + taur*r*AA - (r-n_p)*BB - GG;

    gsl_vector_set (f, 0, y0);
    gsl_vector_set (f, 1, y1);

    return GSL_SUCCESS;
}

struct rparams
{
    double a;
    double b;
};

int
print_state (size_t iter, gsl_multiroot_fsolver * s)
{
printf ("iter = %3u x = % .3f % .3f "
        "f(x) = % .3e % .3e\n",
        iter,
        gsl_vector_get (s->x, 0),
        gsl_vector_get (s->x, 1),
        gsl_vector_get (s->f, 0),
        gsl_vector_get (s->f, 1));
}

int main(){


    // initialize labor supply, pension payments and tax rates
    LL = (2.0+n_p)/(1.0+n_p);
    taup = kappa/((2.0+n_p)*(1.0+n_p));
    tauc = 0.0;
    tauw = 0.0;
    taur = 0.0;

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t i, iter = 0;

    const size_t n = 2;
    struct rparams p = {1.0, 10.0};
    gsl_multiroot_function f = {&eqns_f, n, &p};


    double x_init[2] = {0.1, 0.1};
    gsl_vector *x = gsl_vector_alloc (n);

    gsl_vector_set (x, 0, x_init[0]);
    gsl_vector_set (x, 1, x_init[1]);

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc (T, 2);
    gsl_multiroot_fsolver_set (s, &f, x);

    print_state (iter, s);

    do
        {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);

        print_state (iter, s);

        if (status)   /* check if solver is stuck */
            break;

        status =
            gsl_multiroot_test_residual (s->f, 1e-7);
        }
    while (status == GSL_CONTINUE && iter < 1000);

    printf ("status = %s\n", gsl_strerror (status));

    gsl_multiroot_fsolver_free (s);
    gsl_vector_free (x);
    return 0;
}