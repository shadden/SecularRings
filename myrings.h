#include <stdlib.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include<gsl/gsl_sf.h>
/* orbits.c */
typedef struct {
    double a;
    double e,I,pomega,Omega;
} orbit;

void
OrbitToXYZ(orbit * orb, double u, double XYZ[3]);

void 
OrbitToIntegrationVariables(orbit * orb, double y[4]);

void 
IntegrationVariablesToOrbit(const double y[4], const double a, orbit * orb);

/* averaging.c */
void 
getABC(const double e1, orbit * orb, const double u2,double ABC[4],const double b);

void
lambda_roots(double ABC[4], double *l0, double *l1, double *l2);

void
SinglyAverageGradH_pomega_Omega_ecc_inc(const double e1, orbit * orb, const double u2,double Fav[4]);

void 
DoubleAverageGrad(const double e1, orbit * orb, double Ftot[4], const double eps_abs[4], const double eps_rel[4]);

void 
DoubleAverageForce(const double e1, orbit * orb, double Ftot[4], const double eps_abs, const double eps_rel);

typedef struct {
    size_t neval;
    int code;
    double error;
} grad_e_correction_error;

double 
grad_e_correction(const double e1,orbit * orb, const double rel_tol,const double abs_tol, gsl_integration_workspace * w , size_t wsize, grad_e_correction_error * err);
/* vectors.c */
void printvector(double *X,int N);

void
rotate_x(const double w[3], const double theta, double u[3]);

void 
rotate_z(const double w[3], const double theta, double u[3]);

double
dot(const double x[3], const double y[3]);

double
norm(const double x[3]);
