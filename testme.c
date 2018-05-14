#include <myrings.h>
#include "stdio.h"
typedef struct{
    double a;
    double e1;
    gsl_integration_workspace * w;
    size_t wsize;
    double force_abs_tol,force_rel_tol;
} ode_f_pars;
int ode_f(double t, const double y[], double f[],void * p){
    ode_f_pars * pars;
    pars = (ode_f_pars *) p;
    orbit orb;
    IntegrationVariablesToOrbit(y,pars->a,&orb);
    grad_e_correction_error err;
    DoubleAverageForce(pars->e1,&orb,f,pars->force_abs_tol,pars->force_rel_tol);
    f[0]+=grad_e_correction(pars->e1, &orb, pars->force_abs_tol,pars->force_rel_tol, pars->w, pars->wsize, &err);
    return GSL_SUCCESS;
}

int main(void){
    orbit * o = malloc(sizeof(orbit));
    const double e1 = 0.33;
    double u;
    o->a = 0.25;
    o->e = 0.9;
    o->I = 0.5 * M_PI - 0.3;
    o->Omega = 2.0;
    o->pomega = M_PI;

    double Fav[4];    
    int Npts = 20;
    double pmg;
    const double eps_abs = 1E-5;
    const double eps_rel = 1.0;
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
    double t = 0.0, t1 = 100.0;
    double y[4];
    printf("%.8f\t%.8f\t%.8f\t%.8f\n",o->pomega,o->Omega,o->e,o->I);
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
      printf("%.8f\t%.8f\t%.8f\t%.8f\n",o->pomega,o->Omega,o->e,o->I);

    }
// free your memory!
// 
    free(o);
    gsl_integration_workspace_free (w);
    gsl_odeiv2_driver_free (d);
}

// for(int i=0;i<=Npts;i++){
//     pmg = i * ( 2*M_PI / (double) Npts);
//     o->pomega = pmg;
//     DoubleAverageForce(e1,o,Fav,eps_abs,eps_rel);
//     Fav[0]+=grad_e_correction(e1, o, eps_abs, eps_rel, w, wsize, &err);
//     printvector(Fav,4);
//     printvector(dummy,4);
// }
