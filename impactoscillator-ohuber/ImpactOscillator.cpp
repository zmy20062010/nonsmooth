#include "SiconosKernel.hpp"
#include <math.h>
//#define WITH_FRICTION
//#define DISPLAY_INTER
using namespace std;
#define PI 3.14159265358979323846
int main(int argc, char* argv[])
{
  try{
    // ================= Creation of the model =======================

    // system model parameters definition
	double m = 1.2;
	double m_1 = 2.5;
	double m_2 = 1.8;
	double l = 1.5;
	double Jt = 1.0 / 3 * m * l * l;
	double ke = 2.0e4;
	double Fl = 150.0;
	double f_excite = 200;
	double we = 2 * PI * f_excite;
	double Fa = 50;  // excitation force = Fa * sin(we * t);
	
	// system contact parameters definition
	double eN = 0.2;
	double eT = 0.0;
	double mu = 0.1;
	
	// system initial state parameters definition
	double phi_0 = PI / 4;
	double omega_0 = 0.0;
	double d0 =  l * sin(phi_0);	
	double x1_0 = 0.0;
	double v1_0 = 0.0;
	double x2_0 = 0.0;
	double v2_0 = 0.0;
		
	// system simulation parameters definition
	unsigned int nDof = 3; // degrees of freedom for robot arm
    double t0 = 0;         // initial computation time
    double T = 0.2;       // final computation time
    double h = 1e-5;       // time step : do not decrease, because of strong
	
	// specifying initial conditions
    SP::SiconosVector q0(new SiconosVector(nDof));
	// q0 --> x1, x2, phi;
    SP::SiconosVector v0(new SiconosVector(nDof));
	// v0 --> v1, v2, omega;
    q0->zero();
    v0->zero();
	(*q0)(0) = x1_0;
    (*q0)(1) = x2_0;
    (*q0)(2) = phi_0;
	
	(*v0)(0) = v1_0;
    (*v0)(1) = v2_0;
    (*v0)(2) = omega_0;
	
    // ============================================
    // ============= Dynamical systems ============
    // ============================================
    cout << "====> Model loading ..." << endl << endl;
	
	SP::LagrangianDS oscillator(new LagrangianDS(q0, v0, "ImpactOscillatorPlugin:mass"));
    oscillator->setComputeNNLFunction("ImpactOscillatorPlugin", "NNL");
    oscillator->setComputeJacobianNNLqFunction("ImpactOscillatorPlugin", "jacNNLq");
    oscillator->setComputeJacobianNNLqDotFunction("ImpactOscillatorPlugin", "jacNNLqDot");
	
    oscillator->setComputeFIntFunction("ImpactOscillatorPlugin", "FInt");
    oscillator->setComputeJacobianFIntqFunction("ImpactOscillatorPlugin", "jacFIntq");
    oscillator->setComputeJacobianFIntqDotFunction("ImpactOscillatorPlugin", "jacFIntqDot");
	
	oscillator->setComputeFExtFunction("ImpactOscillatorPlugin", "FExt");
	
    // ============================================
    // ================ Interactions ==============
    // ============================================
    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(eN, eT, mu, 2));
    SP::Relation relation(new LagrangianScleronomousR("ImpactOscillatorPlugin:Gap", "ImpactOscillatorPlugin:Gapjac"));
    SP::Interaction inter(new Interaction(2, nslaw, relation, 1));
	
    // ============================================
    // ================== Model ===================
    // ============================================
    SP::Model impactOscillator(new Model(t0, T));
    impactOscillator->nonSmoothDynamicalSystem()->insertDynamicalSystem(oscillator);
    impactOscillator->nonSmoothDynamicalSystem()->link(inter, oscillator);
	
    // ============================================
    // ============== Simulation ==================
    // ============================================	
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(oscillator, 0.5, 0.0));
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
    SP::OneStepNSProblem osnspb(new FrictionContact(2, SICONOS_FRICTION_2D_ENUM));

    osnspb->numericsSolverOptions()->dparam[0] = 1e-08;
    osnspb->numericsSolverOptions()->iparam[0] = 100;
    osnspb->numericsSolverOptions()->iparam[2] = 1; // random
    SP::TimeStepping s(new TimeStepping(t));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(osnspb, SICONOS_OSNSP_TS_VELOCITY);
    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(200);

    SP::Topology topo = impactOscillator->nonSmoothDynamicalSystem()->topology();
    // ========= End of model definition=============
	
    // ================================= Computation =================================
	
    // ============================================
    // ============== Computation =================
    // ============================================	
    cout << "====> initialization ... It will soon be done!" << endl << endl;
    impactOscillator->initialize(s);
    int N = ceil((T - t0) / h); // Number of time steps

    // ============================================
	// --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    // ============================================		
    unsigned int outputSize = 11;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = oscillator->q();
    SP::SiconosVector v = oscillator->velocity();
    SiconosVector& y = *inter->y(0);
    SiconosVector& lambda = *inter->lambda(1);
	// q :  x1 , x2 , phi;
	// v :  v1 , v2 , omega;
	dataPlot(0, 0) = impactOscillator->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*v)(0);
    dataPlot(0, 5) = (*v)(1);
    dataPlot(0, 6) = (*v)(2);
	
	// --- Time loop ---
    cout << "====> Start computation ... Please wait with patience!" << endl << endl;

    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
	
	while ((s->hasNextEvent()) && (k<= N))
    {
      s->computeOneStep();
      osnspb->setNumericsVerboseMode(0);
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0); 
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);
      dataPlot(k, 7) = y(0);
      dataPlot(k, 8) = y(1);
      dataPlot(k, 9) = lambda(0);
      dataPlot(k, 10) = lambda(1);


      s->processEvents();
      ++show_progress;
      k++;
    }
	cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // ==== Data output to files =====
    cout << "====> Output file writing ... This may take a few while!" << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
	return 1;
  }
  catch (...)
  {
    cout << "Exception caught in ImpactOscillator.cpp" << endl;
	return 1;
  }
  return 0;
}
