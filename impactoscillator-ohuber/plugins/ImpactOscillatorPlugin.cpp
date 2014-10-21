#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  
#include <math.h>
#include <iostream>
#include "SiconosKernel.hpp"
#include "RuntimeException.hpp"

double m = 1.2;
double m_1 = 2.5;
double m_2 = 1.8;
double l = 1.5;
double Jt = 1.0 / 3 * m * l * l;
double k = 2.0e4;
double Fl = 150.0;
double phi_0 = M_PI / 4;
double d0 =  l * sin(phi_0);
double f_excite = 200;
double we = 2 * M_PI * f_excite;
double Fa = 50.0;  // excitation force = Fa * sin(we * t);




// plugins for smooth equations of motion
SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  // columnwise definition of mass matrix
  // q --> x1, x2, phi;
  mass[0] = m + m_1;
  mass[1] = 0.0;
  mass[2] = 1.0 / 2 * m * l * cos(q[2]);

  mass[3] = 0.0;
  mass[4] = m_2;
  mass[5] = 0.0;

  mass[6] = 1.0 / 2 * m * l * cos(q[2]);
  mass[7] = 0.0;
  mass[8] = Jt;
}

SICONOS_EXPORT void NNL(unsigned int sizeOfq, const double *q, const double *v, double *NNL, unsigned int sizeZ, double* z)
{
  // nonlinear inertia terms (negative in h according to LagrangianDS)
  // q --> x1, x2, phi;
  // v --> v1, v2, omega;
  NNL[0] = - 0.5 * m * l * sin(q[2]) * v[2] * v[2];
  NNL[1] = 0.0;
  NNL[2] = 0.0;
}

SICONOS_EXPORT void jacNNLq(unsigned int sizeOfq, const double *q, const double *v, double *NNLq, unsigned int sizeZ, double* z)
{
  // nonlinear inertia terms (negative in h according to LagrangianDS)
  NNLq[0] = 0.0;
  NNLq[1] = 0.0;
  NNLq[2] = 0.0;
  
  NNLq[3] = 0.0;
  NNLq[4] = 0.0;
  NNLq[5] = 0.0;
  
  NNLq[6] = - 0.5 * m * l * cos(q[2]) * v[2] * v[2];
  NNLq[7] = 0.0;
  NNLq[8] = 0.0;
}

SICONOS_EXPORT void jacNNLqDot(unsigned int sizeOfq, const double *q, const double *v, double *NNLqDot, unsigned int sizeZ, double* z)
{
  // nonlinear inertia terms (negative in h according to LagrangianDS)
  NNLqDot[0] = 0.0;
  NNLqDot[1] = 0.0;
  NNLqDot[2] = 0.0;
  
  NNLqDot[3] = 0.0;
  NNLqDot[4] = 0.0;
  NNLqDot[5] = 0.0;
  
  NNLqDot[6] = - m * l * cos(q[2]) * v[2];
  NNLqDot[7] = 0.0;
  NNLqDot[8] = 0.0;
}

SICONOS_EXPORT void FInt(double time, unsigned int sizeOfq, const double *q, const double *v, double *fInt, unsigned int sizeZ, double* z)
{
  // internal forces (negative in h according to LagrangianDS)
  fInt[0] = 0.0;
  fInt[1] = 0.0;
  fInt[2] = k * (l * sin(q[2]) - d0) * l * cos(q[2]);
}

SICONOS_EXPORT void jacFIntq(double time, unsigned int sizeOfq, const double *q, const double *v, double *fIntq, unsigned int sizeZ, double* z)
{
  // internal forces (negative in h according to LagrangianDS)
  fIntq[0] = 0.0;
  fIntq[1] = 0.0;
  fIntq[2] = 0.0;

  fIntq[3] = 0.0;
  fIntq[4] = 0.0;
  fIntq[5] = 0.0;

  fIntq[6] = 0.0;
  fIntq[7] = 0.0;
  fIntq[8] = k * l * (l * cos(q[2]) * cos(q[2]) - l * sin(q[2]) * sin(q[2]) + d0 * sin(q[2]));  
}

SICONOS_EXPORT void jacFIntqDot(double time, unsigned int sizeOfq, const double *q, const double *v, double *fIntqDot, unsigned int sizeZ, double* z)
{
  // internal forces (negative in h according to LagrangianDS)
  fIntqDot[0] = 0.0;
  fIntqDot[1] = 0.0;
  fIntqDot[2] = 0.0;

  fIntqDot[3] = 0.0;
  fIntqDot[4] = 0.0;
  fIntqDot[5] = 0.0;

  fIntqDot[6] = 0.0;
  fIntqDot[7] = 0.0;
  fIntqDot[8] = 0.0;  
}

SICONOS_EXPORT void FExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  // internal forces (negative in h according to LagrangianDS)
  fExt[0] = Fa * sin(we * time);
  fExt[1] = -Fl;
  fExt[2] = 0.0;
}

SICONOS_EXPORT void Gap(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* gap, unsigned int sizeZ, double* z)
{
  gap[0] = d0 - l * sin(q[2]) - q[0]; // normal
  gap[1] = q[1] - l * cos(q[2]); // tangential
}

SICONOS_EXPORT void Gapjac(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* gapjac, unsigned int sizeZ, double* z)
{
  gapjac[0] = -1; // normal
  gapjac[1] = 0.0; // tangential
  gapjac[2] = - l * cos(q[2]); // normal

  gapjac[3] = 0.0; // tangential  
  gapjac[4] = 1.0; // normal
  gapjac[5] = l * sin(q[2]); // tangential  
}
