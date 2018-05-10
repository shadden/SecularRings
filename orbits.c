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
