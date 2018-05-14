#include <myrings.h>
#include "stdio.h"
int main(void){
    orbit * o = malloc(sizeof(orbit));
    const double e1 = 0.33;
    double u;
    // e,I,pomega,Omega: 9.495445e-01	1.217342e+00	3.065128e+00	9.180694e-01
    o->a = 0.25;
    o->e = 9.495445e-01;
    o->I = 1.217342e+00;
    o->Omega = 9.180694e-01;
    o->pomega = 3.065128e+00;
    double pmg;
    const double eps_abs = 1E-5;
    const double eps_rel = 0.1;
    double result;
    size_t wsize = 1000;
    gsl_integration_workspace * w
         = gsl_integration_workspace_alloc (wsize);
    grad_e_correction_error err;
    grad_e_correction(e1,o,eps_abs,0.0 ,w , wsize, &err);
}
