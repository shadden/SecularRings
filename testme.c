#include <myrings.h>
#include "stdio.h"
int main(void){
    orbit * o = malloc(sizeof(orbit));
    const double e1 = 0.2;
    double u;
    o->a = 0.5;
    o->e = 0.1;
    o->I = 0.3;
    o->pomega = 1;
    o->Omega = 2;
    double XYZ[3];
    OrbitToXYZ(o,u,XYZ);
    double ABC[4];
    getABC(e1,o,u,ABC,0.);
    double l0,l1,l2;
    
    //printvector(XYZ,3);
    //printvector(ABC,4);
    lambda_roots(ABC,&l0,&l1,&l2);
    double lvec[3]={l0,l1,l2};
    // printvector(lvec,3);
    double Fav[4],FavC[4];    
    int Npts = 100;
//    for(int i=0;i<=Npts;i++){
//        u = i * ( 2*M_PI / (double) Npts);
//        SinglyAverageGradH_pomega_Omega_ecc_inc(e1, o, u,Fav);
//        printvector(Fav,4);
//    }
    const double eps_abs = 1E-11;
    const double eps_rel = 1E-5;
    DoubleAverageF(e1,o,Fav,eps_rel,eps_abs);
    printvector(Fav,4);
    
    free(o);
}
