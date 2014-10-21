# define pi 3.14159265358979323846
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
double phi_0 = pi / 4;
double v_0 = 0.0;
double omega_0 = 0.0;

// external excitation
long fr = 200;
double we = 2 * pi * fr;
double x_1_0 = 6e-6;