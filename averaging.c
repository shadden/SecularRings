#include "myrings.h"

static double
sign(const double x) { if (x < 0.0) { return -1.0; } else { return 1.0; } }

void getABC(const double e1, orbit * orb, const double u2,double ABC[4],const double b){
    double XYZ[3];
    const double cosphi1 = sqrt(1 - e1*e1);
    OrbitToXYZ(orb,u2,XYZ);
    double A = b*b + 1 + dot(XYZ,XYZ) + 2 * e1 * XYZ[0];
    double Bcos = XYZ[0] + e1;
    double Bsin = cosphi1 * XYZ[1];
    ABC[0] = A;
    ABC[1] = Bsin;
    ABC[2] = Bcos;
    ABC[3] = e1*e1;
}
void
lambda_roots(double ABC[4], double *l0, double *l1, double *l2) {
  int nroots;

  const double A=ABC[0];
  const double Bsine=ABC[1];
  const double Bcose=ABC[2];
  const double C=ABC[3];

  double B2 = Bcose*Bcose + Bsine*Bsine;
  double CmA = C - A;
  double CmA2 = CmA*CmA;
  double CmA3 = CmA2*CmA;
  double disc = B2 - A*C;
  double Q = 1.0/9.0*CmA2 - 1.0/3.0*disc;
  double R = 1.0/27.0*CmA3 - 1.0/6.0*CmA*disc + 0.5*Bsine*Bsine*C;
  double sqrtQ = sqrt(fabs(Q));
  double CmAO3 = CmA/3.0;

  double cosTheta = R / (sqrtQ*sqrtQ*sqrtQ);

  if (cosTheta > 1.0) cosTheta = 1.0;
  if (cosTheta < -1.0) cosTheta = -1.0;

  double theta = acos(cosTheta);
  
  *l0 = -2.0*sqrtQ*cos(theta/3.0 + 2.0*M_PI/3.0) - CmAO3;
  *l1 = -2.0*sqrtQ*cos(theta/3.0 - 2.0*M_PI/3.0) - CmAO3;
  *l2 = -2.0*sqrtQ*cos(theta/3.0) - CmAO3;
}
void
Qmatrix(const double A, const double Bsine, const double Bcose, const double C,
        const double l0, const double l1, const double l2, 
        double Q[3][3]) {
  
  Q[0][0] = sqrt(l0*(l0+C)/((l0-l1)*(l0-l2)));
  Q[0][1] = sqrt(l1*(l1+C)/((l0-l1)*(l1-l2)));
  Q[0][2] = sqrt(fabs(l2)*fabs(l2+C)/((l0-l2)*(l1-l2)));
  Q[1][0] = Bsine*sqrt((l0+C)/(l0*(l0-l1)*(l0-l2)));
  Q[1][1] = Bsine*sqrt((l1+C)/(l1*(l0-l1)*(l1-l2)));
  Q[1][2] = -Bsine*sqrt(fabs(l2+C)/(fabs(l2)*(l0-l2)*(l1-l2)));
  Q[2][0] = Bcose*sqrt(l0/((l0+C)*(l0-l1)*(l0-l2)));
  Q[2][1] = Bcose*sqrt(l1/((l1+C)*(l0-l1)*(l1-l2)));
  Q[2][2] = Bcose*sqrt(fabs(l2)/(fabs(l2+C)*(l0-l2)*(l1-l2)));

  if (fabs((l2+C)/l2) < 1e-3) {
    /* l2 --> C, re-write Q[2][2] to account for loss of accuracy. */
    Q[2][2] = sign(Bcose)*sqrt(fabs(l2)*(l0+C)*(l1+C)/(C*(l0-l2)*(l1-l2)));
  } 

  if (fabs(l2/C) < 1e-3) {
    Q[1][2] = -sign(Bsine)*sqrt((l2+C)*l0*l1/(C*(l0-l1)*(l1-l2)));
  }

  if (fabs(l1/A) < 1e-3) {
    Q[1][1] = sign(Bsine)*sqrt(-l2*(l1+C)*l0/(C*(l0-l1)*(l1-l2)));
  }
  
}

void
UV_from_Q(double Q[3][3], const double e, double U[3], double V[3]) {
  U[0] = Q[0][0]*Q[0][0] - e*Q[0][0]*Q[2][0] + Q[0][2]*Q[0][2] - e*Q[0][2]*Q[2][2];
  U[1] = Q[0][0]*Q[1][0] - e*Q[1][0]*Q[2][0] + Q[0][2]*Q[1][2] - e*Q[1][2]*Q[2][2];
  U[2] = Q[0][0]*Q[2][0] - e*Q[2][0]*Q[2][0] + Q[0][2]*Q[2][2] - e*Q[2][2]*Q[2][2];
  
  V[0] = Q[0][1]*Q[0][1] - e*Q[0][1]*Q[2][1] - Q[0][2]*Q[0][2] + e*Q[0][2]*Q[2][2];
  V[1] = Q[0][1]*Q[1][1] - e*Q[1][1]*Q[2][1] - Q[0][2]*Q[1][2] + e*Q[1][2]*Q[2][2];
  V[2] = Q[0][1]*Q[2][1] - e*Q[2][1]*Q[2][1] - Q[0][2]*Q[2][2] + e*Q[2][2]*Q[2][2];
}


void SinglyAveragedForce(const double e1, orbit * orb, const double u2, double * F0,  double * F1,  double * F2, double * Fav, int Ndim){
    double ABC[4];
    getABC(e1,orb,u2,ABC,0);
    double l0,l1,l2;
    lambda_roots(ABC,&l0, &l1, &l2);
    double Q[3][3];
    Qmatrix(ABC[0],ABC[1],ABC[2],ABC[3],l0,l1,l2,Q);
    double U[3],V[3];
    UV_from_Q(Q,e1,U,V);
    double FU,FV;
    double k2 = (l1-l2) / (l0-l2); 
    double k = sqrt(fabs(k2));
    double Ek = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);
    double Kk = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);
    double prefactor = 2.0 * sqrt(l0-l2) / (l0-l1)/ (l1-l2) / M_PI;
    for( int i=0; i<Ndim; i++){
        FU = U[0] * F0[i] + U[1] * F1[i] + U[2] * F2[i];
        FV = V[0] * F0[i] + V[1] * F1[i] + V[2] * F2[i];
        Fav[i] = prefactor * ((k2*FU+FV)*Ek - (1.0-k2)*FV*Kk );
    }
}
void SinglyAverageGradH_pomega_Omega_ecc_inc(const double e1, orbit * orb, const double u2,double Fav[4]){
    double F0[4];
    double F1[4];
    double F2[4];
    const int Ndim=4;
    const double x0 = (orb->a) * ( cos(u2) - (orb->e) );     
    const double y0 = (orb->a) * ( sin(u2) * sqrt( 1 - (orb->e) * (orb->e)) );
    double xyz[3],_xyz[3];
    double angle;
    const double pomega = orb->pomega;
    const double Omega = orb->Omega;
    const double e = orb->e;
    const double I = orb->I;
    const double omega = pomega - Omega;
    //
    // dr/dpomega
    xyz[0] = x0; xyz[1]=y0; xyz[2]=0;
    angle = omega + 0.5 * M_PI;
    rotate_z(xyz,angle,_xyz);
    rotate_x(_xyz,I,xyz);
    rotate_z(xyz,Omega,_xyz);
    double dxdw = _xyz[0];
    double dydw = _xyz[1];
    double dzdw = _xyz[2];
    //
    // dr/dOmega
    xyz[0] = x0; xyz[1]=y0; xyz[2]=0;
    double M[3][3];
    const double sini = sin(I);
    const double siniBy2Sq = sin(0.5 * I) * sin(0.5 * I);
    const double sin_pomega_2Omega = sin(pomega - 2* Omega);
    const double cos_pomega_2Omega = cos(pomega - 2* Omega);
    //
    M[0][0] =  2 * siniBy2Sq  * sin_pomega_2Omega;
    M[0][1] =  2 * siniBy2Sq  * cos_pomega_2Omega;
    M[0][2] =  cos(Omega) * sini;
    //
    M[1][0] =  2 * siniBy2Sq * cos_pomega_2Omega; 
    M[1][1] = -2 * siniBy2Sq * sin_pomega_2Omega;
    M[1][2] =  2 * sin(Omega) * sini;
    //
    M[2][0] = -1 * sini * cos(omega);
    M[2][1] = sini * sin(omega);
    M[2][2] = 0;
    double dxdOmega = M[0][0] * xyz[0] + M[0][1] * xyz[1] + M[0][2] * xyz[2];
    double dydOmega = M[1][0] * xyz[0] + M[1][1] * xyz[1] + M[1][2] * xyz[2];
    double dzdOmega = M[2][0] * xyz[0] + M[2][1] * xyz[1] + M[2][2] * xyz[2];
    //
    // dr/dI
    xyz[0] = x0; xyz[1]=y0; xyz[2]=0;
    rotate_z(xyz,omega,_xyz);
    rotate_x(_xyz,I+0.5*M_PI,xyz);
    xyz[0]=0;
    rotate_z(xyz,orb->Omega,_xyz);
    double dxdI = _xyz[0];
    double dydI = _xyz[1];
    double dzdI = _xyz[2];
    //
    // dr/de
    xyz[0] = -1 * orb->a ; xyz[1]= -1 * orb->a * e * sin(u2) / sqrt(1-e*e) ; xyz[2]=0;
    rotate_z(xyz,omega,_xyz);
    rotate_x(_xyz,I,xyz);
    rotate_z(xyz,Omega,_xyz);
    double dxde = _xyz[0];
    double dyde = _xyz[1];
    double dzde = _xyz[2];
    /*******************************/
    // F0 = r.grad(r) - e1*grad(x)
    // F1 = -sqrt(1-e1*e1)*grad(y)
    // F2 = -grad(x) 
    //
    // Note r.(dr/dpomega) = r.(dr/dOmega) = r.(dr/dI) =0
    /*******************************/
    OrbitToXYZ(orb,u2,xyz);
    const double x = xyz[0];
    const double y = xyz[1];
    const double z = xyz[2];
// F0
    F0[0] =  e1 * dxdw ; 
    F0[1] =  e1 * dxdOmega ; 
    F0[2] =  x*dxde + y*dyde + z*dzde + e1 * dxde ; 
    F0[3] =  e1 * dxdI ; 
    
// F1
    double cosphi1 = -1*sqrt(1-e1*e1);
    F1[0] = cosphi1 * dydw;
    F1[1] = cosphi1 * dydOmega;
    F1[2] = cosphi1 * dyde;
    F1[3] = cosphi1 * dydI;
// F2
    F2[0] = -1 * dxdw;
    F2[1] = -1 * dxdOmega;
    F2[2] = -1 * dxde;
    F2[3] = -1 * dxdI;
    
    SinglyAveragedForce(e1,orb,u2,F0,F1,F2,Fav,Ndim);
    
}
void DoubleAverageGrad(const double e1, orbit * orb, double Ftot[4], const double eps_abs[4], const double eps_rel[4]){
    double F0[4],FC[4];
    double Ftot_old[4];
    double abs_err[4], rel_err[4];
    double max_abs_err, max_rel_err;
    int N,i,j;
    const int NMax = 1048576;
    const double e = orb->e;

    memset(Ftot, 0, 4*sizeof(double));
    memset(Ftot_old, 0, 4*sizeof(double));
    memset(F0, 0, 4*sizeof(double));
    memset(FC, 0, 4*sizeof(double));
    /* Compute the u = 0 term. */
    double F[4];
    bool tol = false;
    bool tolj;
    SinglyAverageGradH_pomega_Omega_ecc_inc(e1,orb,0.0,F);
    for (i = 0; i < 4; i++) {
      F0[i] += F[i];
      FC[i] += F[i];
    }
    N = 1;
    do {
      memcpy(Ftot_old,Ftot,4*sizeof(double));
      double h;
      N *= 2;
      h = 1.0 / N;
      for (i = 0; i < 4; i++) {
        F0[i] /= 2.0;
        FC[i] /= 2.0;
      }
      for (i = 1; i < N; i += 2) {
        double u = 2.0*M_PI*((double) i)/((double) N);
        double c = cos(u);
        SinglyAverageGradH_pomega_Omega_ecc_inc(e1,orb,u,F);
        for(j=0; j<4; j++){
            F0[j] += h*F[j];
            FC[j] += c*h*F[j];
        }
     }
     tol = true;
     tolj = true;
     for(j=0; j<4; j++){
        Ftot[j] = F0[j] - e * FC[j];
        abs_err[j] = fabs(Ftot[j]-Ftot_old[j]);
        rel_err[j] = fabs(abs_err[j] / Ftot[j]);
        tolj = ((rel_err[j] < eps_rel[j]) || (abs_err[j] < eps_abs[j]));
        tol = tol && tolj;
     }
    
    } while(( N<16 || !(tol) ) && (N<NMax)); //((N < 16 || err > epsacc) && (N < NMax) && (derr > 0.1));

}
void 
DoubleAverageForce(const double e1, orbit * orb, double Force[4], const double eps_abs, const double eps_rel){
    const double a = orb->a;
    const double I = orb->I;
    const double e = orb->e;
    // d(pomega)/dt = -[ sqrt(1-e^2) / ( sqrt(alpha) * e ) ] * dH/de
    // d(Omega)/dt =  -1 / ( sqrt(alpha) * e )]  dH/dI
    // dG/dt =  -dH/dpomega
    // dZ/dt =  dH/dOmega
    const double pmgdot_factor = -1.0 * sqrt(1.0-e*e) / sqrt(a) / e ;
    const double Omgdot_factor = -1.0 / sqrt(1.0-e*e) / sqrt(a) / sin(I) ; 
    const double Gdot_factor = -1.0; 
    const double Zdot_factor = 1.0; 
    // 0 - dH/dpomega
    // 1 - dH/dOmega
    // 2 - dH/de
    // 3 - dH/dI
    double Grad[4],eps_rel_arr[4],eps_abs_arr[4];
    
    for(int i=0;i<4;i++){eps_rel_arr[i]=eps_rel;eps_abs_arr[i]=eps_abs;}
    eps_abs_arr[0]*=fabs(pmgdot_factor);
    eps_abs_arr[1]*=fabs(Omgdot_factor);

    DoubleAverageGrad(e1, orb, Grad, eps_abs_arr, eps_rel_arr);
    Force[0] = pmgdot_factor * (Grad[2]);
    Force[1] = Omgdot_factor * Grad[3];
    Force[2] = Gdot_factor * Grad[0];
    Force[3] = Zdot_factor * Grad[1];
    
    
}

typedef struct  {
 orbit * o;
 double e1,u2;
} inner_u1_integrand_params;

typedef struct  {
 orbit * o;
 double e1,rel_tol,abs_tol,error;
 size_t neval, wsize;
 int code;
 gsl_integration_workspace * w;
} outer_u2_integrand_params;

double 
inner_u1_integrand(double u1, void * p){
    inner_u1_integrand_params * params;
    params = (inner_u1_integrand_params *) p;
    const double e1 = params->e1;
    double XYZ[3];
    OrbitToXYZ(params->o,params->u2,XYZ);
    const double dx = cos(u1) - e1 - XYZ[0]; 
    const double dy = sqrt(1.0-e1*e1) * sin(u1) - XYZ[1];
    const double dz = XYZ[2];
    const double Delta = sqrt(dx*dx+dy*dy+dz*dz);
    return  (1.0 - e1 * cos(u1) )  / Delta / 2.0 / M_PI;
}
double 
outer_u2_integrand(double u2, void * p){
    outer_u2_integrand_params * outer_params;
    outer_params = (outer_u2_integrand_params *) p;
    inner_u1_integrand_params inner_params;
    inner_params.o = outer_params->o;
    inner_params.e1= outer_params->e1;
    inner_params.u2 = u2;
    gsl_function F;
    F.function = &inner_u1_integrand;
    F.params = (void *)(&inner_params);
    double result;
    // outer_params->code=gsl_integration_qng(&F,0.0,2.0*M_PI,(outer_params->abs_tol),(outer_params->rel_tol),&result,&(outer_params->error),&(outer_params->neval));
    outer_params->code=gsl_integration_qags(&F,0.0,2.0*M_PI,(outer_params->abs_tol),(outer_params->rel_tol),(outer_params->wsize),(outer_params->w),&result,&(outer_params->error));
    return result * cos(u2) / 2.0 / M_PI;
}

double grad_e_correction(const double e1,orbit * orb, const double abs_tol, const double rel_tol, gsl_integration_workspace * w , size_t wsize , grad_e_correction_error * err){
    double result;
    outer_u2_integrand_params pars;
    pars.o= orb;
    pars.e1=e1;
    pars.abs_tol=abs_tol;
    pars.rel_tol=rel_tol;
    pars.w = w;
    pars.wsize=wsize;
    gsl_function F;
    F.function = &outer_u2_integrand;
    F.params = (void *) &pars;
    //err->code = gsl_integration_qng(&F,0.0,2*M_PI, abs_tol, rel_tol, &result,&(err->error),&(err->neval));
    err->code=gsl_integration_qags(&F,0.0,2.0*M_PI, abs_tol, rel_tol,wsize,w,&result,&(err->error));
    const double e2 = orb->e;
    const double a  = orb->a;
    const double pmgfactor = -1 * sqrt(1-e2*e2) / e2 / sqrt(a);
    return pmgfactor * result;
}

