#include <myrings.h>
#include "stdio.h"
typedef struct{
    double a;
    double e1;
    gsl_integration_workspace * w;
    size_t wsize;
    double force_abs_tol,force_rel_tol;
} ode_f_pars;
void ode_f(double t, const double y[], double f[],void * p){
    ode_f_pars * pars;
    pars = (ode_f_pars *) p;
    orbit orb;
    IntegrationVariablesToOrbit(y,pars->a,&orb);
    grad_e_correction_error err;
    DoubleAverageForce(pars->e1,orb,f,pars->force_abs_tol,pars->force_rel_tol);
    f[0]+=grad_e_correction(pars->e1, orb, pars->force_abs_tol,pars->force_rel_tol, w, wsize, &err);
}

int main(void){
    orbit * o = malloc(sizeof(orbit));
    const double e1 = 0.03;
    double u;
    o->a = 0.5;
    o->e = 0.05;
    o->I = 0.1;
    
    o->Omega = 2;

    double Fav[4];    
    double dummy[4];    
    dummy[0]=dummy[1]=dummy[2]=dummy[3]=0.0;
    int Npts = 20;
    double pmg;
    const double eps_abs = 1E-9;
    const double eps_rel = 1e-2;
    double result;
    size_t wsize = 1000;
    gsl_integration_workspace * w
         = gsl_integration_workspace_alloc (wsize);
    grad_e_correction_error err;

    for(int i=0;i<=Npts;i++){
        pmg = i * ( 2*M_PI / (double) Npts);
        o->pomega = pmg;
        DoubleAverageForce(e1,o,Fav,eps_abs,eps_rel);
        Fav[0]+=grad_e_correction(e1, o, eps_abs, eps_rel, w, wsize, &err);
        printvector(Fav,4);
        printvector(dummy,4);
    }
// ODE integration!!
    void * jac; // null pointer for jacobian as placeholder
    ode_f_pars pars;
    pars.w = w;
    pars.wsize=wsize;
    pars.force_abs_tol=1.0E-7;
    pars.force_rel_tol=1.0E-2;
    pars.e1 = e1;
    pars.a = o.a
    gsl_odeiv2_system sys {ode_f,jac,pars};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk45,
                                  1e-6, 1e-6, 0.0);
// free your memory!
// 
    free(o);
    gsl_integration_workspace_free (w);
    gsl_odeiv2_driver_free (d);
}

