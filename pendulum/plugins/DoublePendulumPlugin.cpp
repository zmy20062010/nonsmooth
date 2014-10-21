#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <stdio.h>
#include <math.h>

double g = 10.0;
double m1 = 1.0;
double m2 = 1.0 ;
double l1 = 1.0 ;
double l2 = 1.0 ;

SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  mass[0] = (m1 + m2) * l1;
  mass[1] = m2 * l1 * cos(q[0] - q[1]);
  mass[2] = m2 * l2 * cos(q[0] - q[1]);
  mass[3] = m2 * l2;
}

SICONOS_EXPORT void NNL(unsigned int sizeOfq, const double *q, const double *v, double *NNL, unsigned int sizeZ, double* z)
{
  NNL[0] =    m2 * l2 * v[1] * v[1] * sin(q[0] - q[1]);
  NNL[1] =  - m2 * l1 * v[0] * v[0] * sin(q[0] - q[1]);
}

SICONOS_EXPORT void jacobianNNLq(unsigned int sizeOfq, const double *q, const double *v, double *jacob, unsigned int sizeZ, double* z)
{
  jacob[0] =  m2 * l2 * v[1] * v[1] * cos(q[0] - q[1]);
  jacob[1] =  - m2 * l1 * v[0] * v[0] * cos(q[0] - q[1]);
  jacob[2] =  - m2 * l2 * v[1] * v[1] * cos(q[0] - q[1]);
  jacob[3] =  m2 * l1 * v[0] * v[0] * cos(q[0] - q[1]);
}

SICONOS_EXPORT void jacobianVNNL(unsigned int sizeOfq, const double *q, const  double *v, double *jacob, unsigned int sizeZ, double* z)
{
  jacob[0] =  0.0;
  jacob[1] =  - 2.0 * m2 * l1 * v[0] * sin(q[0] - q[1]);
  jacob[2] =  2.0 * m2 * l2 * v[1] * sin(q[0] - q[1]);
  jacob[3] =  0.0;
}

SICONOS_EXPORT void FInt(double time, unsigned int sizeOfq, const double *q, const double *v, double *fInt, unsigned int sizeZ, double* z)
{
  fInt[0] = (m1 + m2) * g * sin(q[0]);
  fInt[1] = m2 * g * sin(q[1]);
}
SICONOS_EXPORT void jacobianFIntq(double time, unsigned int sizeOfq, const double *q, const double *v, double *jacob, unsigned int sizeZ, double* z)
{
  jacob[0] = (m1 + m2) * g * cos(q[0]);
  jacob[1] = 0.0;
  jacob[2] = 0.0
  jacob[3] = m2 * g * cos(q[1]);
}

SICONOS_EXPORT void jacobianVFInt(double time, unsigned int sizeOfq, const double *q, const double *v, double *jacob, unsigned int sizeZ, double* z)
{
  jacob[0] = 0.0;
  jacob[1] = 0.0;
  jacob[2] = 0.0;
  jacob[3] = 0.0;
}

SICONOS_EXPORT void h0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  y[0] = l1 * sin(q[0]);
}

SICONOS_EXPORT void G0(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  G[0] = l1 * cos(q[0]);
  G[1] = 0.0;
}

SICONOS_EXPORT void h1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* y, unsigned int sizeZ, double* z)
{
  y[0] = l1 * sin(q[0]) + l2 * sin(q[1]);
}

SICONOS_EXPORT void G1(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* G, unsigned int sizeZ, double* z)
{
  G[0] = l1 * cos(q[0]);
  G[1] = l2 * cos(q[1]);
}

