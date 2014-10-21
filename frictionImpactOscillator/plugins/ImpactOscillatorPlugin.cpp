#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <stdio.h>
#include <math.h>
#include "frictionImpactOscillator.h"


SICONOS_EXPORT void FInt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fInt, unsigned int sizeZ, double* z)
{
  fInt[0] = 0.0;
  fInt[1] = (k * l * (sin(phi_0) - sin(q[1])) - 1. / 2 * m * we * we * x_1_0 * cos(we * time)) * l * cos(q[1]);
}

SICONOS_EXPORT void FExt(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *fExt, unsigned int sizeZ, double* z)
{
  fExt[0] = - Fl;
  fExt[1] = 0.0;
}


SICONOS_EXPORT void h(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *y, unsigned int sizeZ, double* z)
{
  y[0] = l * sin(phi_0) - l * sin(q[1]);
  y[1] = q[0] - l * cos(q[1]);
}

SICONOS_EXPORT void jach((double time, unsigned int sizeOfq, const double *q, const double *velocity, double *dy, unsigned int sizeZ, double* z)
{
   dy(0,0) = 0.0;
   dy(0,1) = - l * cos(q[1]) ;
   dy(1,0) = 1.0;
   dy(1,1) = l * sin(q[1]);
}

SICONOS_EXPORT void doth(double time, unsigned int sizeOfq, const double *q, const double *velocity, double *ht, unsigned int sizeZ, double* z)
{
  ht[0] = we * x_1_0 * sin(we * time);
  ht[1] = 0.0;
}


