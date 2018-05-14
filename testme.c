#include <myrings.h>
#include "stdio.h"
// function parameters for use with 'ode_f'
typedef struct{
    double a;
    double e1;
    gsl_integration_workspace * w;
    size_t wsize;
    double force_abs_tol,force_rel_tol;
} ode_f_pars;
//
// function giving dydt in form used by gsl's ode solvers 
int ode_f(double t, const double y[], double f[],void * p){
    ode_f_pars * pars;
    pars = (ode_f_pars *) p;
    orbit orb;
    IntegrationVariablesToOrbit(y,pars->a,&orb);
    grad_e_correction_error err;
    DoubleAverageForce(pars->e1,&orb,f,pars->force_abs_tol,pars->force_rel_tol);
    // 
    const double abs_err = fmin( fabs(pars->force_rel_tol * f[0]) , pars->force_abs_tol );
    f[0]+=grad_e_correction(pars->e1, &orb,abs_err, 0.0 , pars->w, pars->wsize, &err);
    return GSL_SUCCESS;
}

int main(void){
    orbit * o = malloc(sizeof(orbit));
    const double e1 = 0.33;
    double u;
    o->a = 0.15;
    o->e = 0.05;
    o->I = 0.1;
    o->Omega = 2.0;
    o->pomega = 0.0;

    double Fav[4];    
    int Npts = 20;
    double pmg;
    const double eps_abs = 1E-12;
    const double eps_rel = 1E-5;
    double result;
    size_t wsize = 1000;
    gsl_integration_workspace * w
         = gsl_integration_workspace_alloc (wsize);
    grad_e_correction_error err;
// ODE integration!!
    void * jac; // null pointer for jacobian as placeholder
    ode_f_pars pars;
    pars.w = w;
    pars.wsize=wsize;
    pars.force_abs_tol=eps_abs;
    pars.force_rel_tol=eps_rel;
    pars.e1 = e1;
    pars.a = o->a;
    size_t dim = 4;
    gsl_odeiv2_system sys = {ode_f,jac,dim,(void *) &pars};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
                                  1e-6, 1e-6, 0.0);
    int i;
    double t = 0.0, t1 = 200.0;
    double y[4];
    printf("%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",0.0,o->pomega,o->Omega,o->e,o->I);
    OrbitToIntegrationVariables(o,y);
    for (i = 1; i <= 100; i++)
    {
      double ti = i * t1 / 100.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
      if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }

      IntegrationVariablesToOrbit(y,o->a,o);
      printf("%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n",ti,o->pomega,o->Omega,o->e,o->I);


    }
// free your memory!
// 
    free(o);
    gsl_integration_workspace_free (w);
    gsl_odeiv2_driver_free (d);
}
