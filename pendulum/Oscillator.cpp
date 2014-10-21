/* Siconos-sample , Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */


/*!\file
  C++ input file, MoreauJeanOSI-Time-Stepping version
  T. Schindler, V. Acary

  Slider-crank simulation with a MoreauJeanOSI-Time-Stepping scheme

  see Flores/Leine/Glocker : Modeling and analysis of planar rigid multibody systems with
  translational clearance joints based on the non-smooth dynamics approach
  */

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

int main(int argc, char* argv[])
{
  try
  {
    // ================= Creation of the model =======================

    // parameters according to Table 1
    unsigned int nDof = 3; // degrees of freedom for robot arm
    double t0 = 0;         // initial computation time
    double T = 0.2;       // final computation time
    double h = 1e-5;       // time step : do not decrease, because of strong penetrations

    // initial conditions
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    q0->zero();
    (*q0)(0) = y_0;
    (*q0)(1) = phi_0;	
    v0->zero();
    (*v0)(0) = v_0;
    (*v0)(1) = omega_0;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ..." << endl << endl;

    SP::LagrangianDS slider(new LagrangianDS(q0, v0, "OscillatorPlugin:mass"));
    slider->setComputeNNLFunction("OscillatorPlugin", "NNL");
    slider->setComputeJacobianNNLqFunction("OscillatorPlugin", "jacobianNNLq");
    slider->setComputeJacobianNNLqDotFunction("OscillatorPlugin", "jacobianNNLqDot");
    slider->setComputeFIntFunction("OscillatorPlugin", "FInt");
    slider->setComputeJacobianFIntqFunction("OscillatorPlugin", "jacobianFIntq");
    slider->setComputeJacobianFIntqDotFunction("OscillatorPlugin", "jacobianFIntqDot");

    // -------------------
    // --- Interactions---
    // -------------------
    // -- corner 1 --
#ifdef WITH_FRICTION
    SP::NonSmoothLaw nslaw1(new NewtonImpactFrictionNSL(eN1, eT1, mu1, 2));
    SP::Relation relation1(new LagrangianScleronomousR("SliderCrankPlugin:g1", "SliderCrankPlugin:W1"));
    SP::Interaction inter1(new Interaction(2, nslaw1, relation1, 1));

    // -- corner 2 --
    SP::NonSmoothLaw nslaw2(new NewtonImpactFrictionNSL(eN2, eT2, mu2, 2));
    SP::Relation relation2(new LagrangianScleronomousR("SliderCrankPlugin:g2", "SliderCrankPlugin:W2"));
    SP::Interaction inter2(new Interaction(2, nslaw2, relation2, 2));

    // -- corner 3 --
    SP::NonSmoothLaw nslaw3(new NewtonImpactFrictionNSL(eN3, eT3, mu3, 2));
    SP::Relation relation3(new LagrangianScleronomousR("SliderCrankPlugin:g3", "SliderCrankPlugin:W3"));
    SP::Interaction inter3(new Interaction(2, nslaw3, relation3, 3));

    // -- corner 4 --
    SP::NonSmoothLaw nslaw4(new NewtonImpactFrictionNSL(eN4, eT4, mu4, 2));
    SP::Relation relation4(new LagrangianScleronomousR("SliderCrankPlugin:g4", "SliderCrankPlugin:W4"));
    SP::Interaction inter4(new Interaction(2, nslaw4, relation4, 4));
#else
    // -- corner 1 --
    SP::NonSmoothLaw nslaw1(new NewtonImpactNSL(eN1));
    SP::Relation relation1(new LagrangianScleronomousR("SliderCrankPlugin:g1", "SliderCrankPlugin:W1"));
    SP::Interaction inter1(new Interaction(1, nslaw1, relation1, 1));

    // -- corner 2 --
    SP::NonSmoothLaw nslaw2(new NewtonImpactNSL(eN2));
    SP::Relation relation2(new LagrangianScleronomousR("SliderCrankPlugin:g2", "SliderCrankPlugin:W2"));
    SP::Interaction inter2(new Interaction(1, nslaw2, relation2, 2));

    // -- corner 3 --
    SP::NonSmoothLaw nslaw3(new NewtonImpactNSL(eN3));
    SP::Relation relation3(new LagrangianScleronomousR("SliderCrankPlugin:g3", "SliderCrankPlugin:W3"));
    SP::Interaction inter3(new Interaction(1, nslaw3, relation3, 3));

    // -- corner 4 --
    SP::NonSmoothLaw nslaw4(new NewtonImpactNSL(eN4));
    SP::Relation relation4(new LagrangianScleronomousR("SliderCrankPlugin:g4", "SliderCrankPlugin:W4"));
    SP::Interaction inter4(new Interaction(1, nslaw4, relation4, 4));
#endif

    // -------------
    // --- Model ---
    // -------------
    SP::Model sliderWithClearance(new Model(t0, T));
    sliderWithClearance->nonSmoothDynamicalSystem()->insertDynamicalSystem(slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter1, slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter2, slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter3, slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter4, slider);

    // ----------------
    // --- Simulation ---
    // ----------------
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(slider, 0.5, 0.0));
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
#ifdef WITH_FRICTION
    SP::OneStepNSProblem impact(new FrictionContact(2, SICONOS_FRICTION_2D_ENUM));
#else
    SP::OneStepNSProblem impact(new LCP(SICONOS_LCP_LEMKE));
#endif
    impact->numericsSolverOptions()->dparam[0] = 1e-08;
    impact->numericsSolverOptions()->iparam[0] = 100;
    impact->numericsSolverOptions()->iparam[2] = 1; // random
    SP::TimeStepping s(new TimeStepping(t));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_TS_VELOCITY);
    s->setNewtonTolerance(1e-10);
    s->setNewtonMaxIteration(200);

    SP::Topology topo = sliderWithClearance->nonSmoothDynamicalSystem()->topology();
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    cout << "====> Initialisation ..." << endl << endl;
    sliderWithClearance->initialize(s);
    int N = ceil((T - t0) / h) + 1; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 27;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = slider->q();
    SP::SiconosVector v = slider->velocity();

    dataPlot(0, 0) = sliderWithClearance->t0();
    dataPlot(0, 1) = (*q)(0) / (2.*M_PI); // crank revolution
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*v)(0);
    dataPlot(0, 5) = (*v)(1);
    dataPlot(0, 6) = (*v)(2);
    dataPlot(0, 7) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 1 (normalized)
    dataPlot(0, 8) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 2 (normalized)
    dataPlot(0, 9) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (-c); // y corner 3 (normalized)
    dataPlot(0, 10) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (-c); // y corner 4 (normalized)
    dataPlot(0, 11) = (l1 * cos((*q)(0)) + l2 * cos((*q)(1)) - l2) / l1; // x slider (normalized)
    dataPlot(0, 12) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1))) / c; // y slider (normalized
    dataPlot(0, 13) = (*inter1->y(0))(0) ; // g1
    dataPlot(0, 14) = (*inter2->y(0))(0) ; // g2
    dataPlot(0, 15) = (*inter3->y(0))(0) ; // g3
    dataPlot(0, 16) = (*inter4->y(0))(0) ; // g4
    dataPlot(0, 17) = (*inter1->y(1))(0) ; // dot g1
    dataPlot(0, 18) = (*inter2->y(1))(0) ; // dot g2
    dataPlot(0, 19) = (*inter3->y(1))(0) ; // dot g3
    dataPlot(0, 20) = (*inter4->y(1))(0) ; // dot g4
    dataPlot(0, 21) = (*inter1->lambda(1))(0) ; // lambda1
    dataPlot(0, 22) = (*inter2->lambda(1))(0) ; // lambda1
    dataPlot(0, 23) = (*inter3->lambda(1))(0) ; // lambda3
    dataPlot(0, 24) = (*inter4->lambda(1))(0) ; // lambda4
    dataPlot(0, 25) = 0;
    dataPlot(0, 26) = 0;

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;

    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    SP::InteractionsGraph indexSet1 = topo->indexSet(1);

    while ((s->hasNextEvent()) && (k<= 3000))
//    while ((s->hasNextEvent()))
    {

      std::cout <<"=====================================================" <<std::endl;
      std::cout <<"=====================================================" <<std::endl;
      std::cout <<"=====================================================" <<std::endl;
      std::cout <<"Iteration k = " << k <<std::endl;
      std::cout <<"s->nextTime() = " <<s->nextTime()  <<std::endl;
      std::cout <<"=====================================================" <<std::endl;



      //std::cout << "=============== Step k ="<< k<< std::endl;
      s->advanceToEvent();
      impact->setNumericsVerboseMode(0);
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0) / (2.*M_PI); // crank revolution
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);
      dataPlot(k, 7) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 1 (normalized)
      dataPlot(k, 8) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 2 (normalized)
      dataPlot(k, 9) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (c); // y corner 3 (normalized)
      dataPlot(k, 10) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (c); // y corner 4 (normalized)
      dataPlot(k, 11) = (l1 * cos((*q)(0)) + l2 * cos((*q)(1)) - l2) / l1; // x slider (normalized)
      dataPlot(k, 12) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1))) / c; // y slider (normalized)
      dataPlot(k, 13) = (*inter1->y(0))(0) ; // g1
      dataPlot(k, 14) = (*inter2->y(0))(0) ; // g2
      dataPlot(k, 15) = (*inter3->y(0))(0) ; // g3
      dataPlot(k, 16) = (*inter4->y(0))(0) ; // g4
      dataPlot(k, 17) = (*inter1->y(1))(0) ; // dot g1
      dataPlot(k, 18) = (*inter2->y(1))(0) ; // dot g2
      dataPlot(k, 19) = (*inter3->y(1))(0) ; // dot g3
      dataPlot(k, 20) = (*inter4->y(1))(0) ; // dot g4
      dataPlot(k, 21) = (*inter1->lambda(1))(0) ; // lambda1
      dataPlot(k, 22) = (*inter2->lambda(1))(0) ; // lambda1
      dataPlot(k, 23) = (*inter3->lambda(1))(0) ; // lambda3
      dataPlot(k, 24) = (*inter4->lambda(1))(0) ; // lambda4
      dataPlot(k, 25) = s->getNewtonNbSteps();
      dataPlot(k, 26) = indexSet1->size();

      if (indexSet1->size() > 5)
      {
        impact->display();
      }
      //      if (s->nextTime() > 0.035 and (*inter1->lambda(1))(0) >0.0)
#ifdef DISPLAY_INTER
        std::cout << "=============== Step k =" << k << std::endl;
        std::cout << "Time " << s->nextTime() << std::endl;

        impact->display();
        std::cout << " (*inter1->lambda(1))(0) " << (*inter1->lambda(1))(0) << std:: endl;
        std::cout << " (*inter2->lambda(1))(0) " << (*inter2->lambda(1))(0) << std:: endl;
        std::cout << " (*inter3->lambda(1))(0) " << (*inter3->lambda(1))(0) << std:: endl;
        std::cout << " (*inter4->lambda(1))(0) " << (*inter4->lambda(1))(0) << std:: endl;
#endif

      s->processEvents();
      ++show_progress;
      k++;
    }

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in SliderCrankD1MinusLinearOSI.cpp" << endl;
  }
}