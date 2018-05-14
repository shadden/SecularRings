#include "myrings.h"
#include <math.h>
void OrbitToXYZ(orbit * orb, double u, double XYZ[3]){
    const double sinphi = orb->e;
    const double cosphi = sqrt( 1 - sinphi * sinphi);
    const double x = orb->a * (cos(u) - sinphi) ;
    const double y = orb->a * cosphi * sin(u) ;
    double xyz[3] = {x,y,0};
    double _xyz[3] ;
    const double omega = orb->pomega - orb->Omega;
    rotate_z(xyz,omega,_xyz);
    rotate_x(_xyz,orb->I,xyz);
    rotate_z(xyz,orb->Omega,XYZ);
}
void OrbitToIntegrationVariables(orbit * orb, double y[4]){
    const double a = orb->a;
    const double e = orb->e;
    const double I = orb->I;
    const double G = sqrt(a * (1 - e*e) );
    // y = ( pomega, -Omega, G , Z )
    y[0] = orb->pomega;
    y[1] = -1 * orb->Omega;
    y[2] = G;
    y[3] = G * (1 - cos(I));
}
void IntegrationVariablesToOrbit(const double y[4], const double a, orbit * orb){
    // y = ( pomega, -Omega, G , Z )
    orb->pomega = y[0];
    orb->Omega = -1 * y[1];
    orb->a = a;
    const double G = y[2];
    const double Z = y[3];
    orb->e = sqrt( 1.0 - G*G / a);
    orb->I = acos( 1.0 - Z / G);
}
