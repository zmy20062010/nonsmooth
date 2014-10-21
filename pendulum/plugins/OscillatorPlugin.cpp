#include "SiconosKernel.hpp"
//#define WITH_FRICTION
//#define DISPLAY_INTER
using namespace std;


# define PI 3.14159265358979323846
// geometry
double l = 1.5e-3;

// inertias
double m_2 = 2e-4;
double m = 1.5e-4;
double J = 1.0 / 3 * m * l * l;

// contact parameters
double eps_N = 0.5;
double eps_T = 0;
double mu = 0.3;

// forces
double Fl = 1.0e-6;

// initial conditions
double y_0 = 0.0;
double phi_0 = PI / 4;
double v_0 = 0.0;
double omega_0 = 0.0;

// external excitation
long fr = 200;
double we = 2 * pi * fr;
double x_1_0 = 6e-6;

// q[0] = x2; q[1] = phi;
SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  mass[0] = m_2;
  mass[1] = 0.0;
  mass[2] = 0.0;
  mass[3] = J;
}
SICONOS_EXPORT void NNL(unsigned int sizeOfq, const double *q, const double *v, double *NNL, unsigned int sizeZ, double* z)
{
  NNL[0] =    0.0;
  NNL[1] =    k * l * l * (sin(phi_0) - sin(q[1])) * cos(q[1]);
}

SICONOS_EXPORT void jacobianNNLq(unsigned int sizeOfq, const double *q, const double *v, double *NNLq, unsigned int sizeZ, double* z)
{
  NNLq[0] =  0.0;
  NNLq[1] =  0.0;
  NNLq[2] =  0.0;
  NNLq[3] =  k * l * l * (- sin(phi_0) * cos(q[1]) + sin(q[1]) * sin(q[1]) - cos(q[1]) * cos(q[1]));
}

SICONOS_EXPORT void jacobianNNLqDot(unsigned int sizeOfq, const double *q, const  double *v, double *NNLqDot, unsigned int sizeZ, double* z)
{
  NNLqDot[0] =  0.0;
  NNLqDot[1] =  0.0;
  NNLqDot[2] =  0.0;
  NNLqDot[3] =  0.0;
}

SICONOS_EXPORT void FInt(double time, unsigned int sizeOfq, const double *q, const double *v, double *fInt, unsigned int sizeZ, double* z)
{
  fInt[0] = 0.0;
  fInt[1] = - 0.5 * m * l * x_1_0 * we * we * cos(we * time) * cos(q[1]);
}

SICONOS_EXPORT void jacobianFIntq(double time, unsigned int sizeOfq, const double *q, const double *v, double *FIntq, unsigned int sizeZ, double* z)
{
  FIntq[0] = 0.0;
  FIntq[1] = 0.0;
  FIntq[2] = 0.0
  FIntq[3] = 0.5 * m * l * x_1_0 * we * we * cos(we * time) * sin(q[1]);
}

SICONOS_EXPORT void jacobianFIntqDot(double time, unsigned int sizeOfq, const double *q, const double *v, double *FIntqDot, unsigned int sizeZ, double* z)
{
  FIntqDot[0] = 0.0;
  FIntqDot[1] = 0.0;
  FIntqDot[2] = 0.0;
  FIntqDot[3] = 0.0;
}

SICONOS_EXPORT void FExt(double time, unsigned int sizeOfq, const double *q, const double *v, double *fExt, unsigned int sizeZ, double* z)
{
  fExt[0] = - Fl;
  fExt[1] = 0.0;
}

SICONOS_EXPORT void Gap(double time, unsigned int sizeOfq, const double *q, const double *v, double *gap, unsigned int sizeZ, double* z)
{
  gap[0] = l * sin(phi_0) - x_1_0 + x_1_0 * cos(we * time) - l * sin(q[1]);
  gap[1] = q[0] - l * cos(q[1]);
}

SICONOS_EXPORT void Gapq(double time, unsigned int sizeOfq, const double *q, const double *v, double *gapq, unsigned int sizeZ, double* z)
{
  gapq[0] = 0.0;
  gapq[1] = 1.0;
  gapq[2] = - l * cos(q[1]);
  gapq[3] = l * sin(q[1]);  
}

SICONOS_EXPORT void Gapt(double time, unsigned int sizeOfq, const double *q, const double *v, double *gapt, unsigned int sizeZ, double* z)
{
  gapt[0] = - x_1_0 * we * sin(we * time);
  gapt[1] = 0.0;
}





















