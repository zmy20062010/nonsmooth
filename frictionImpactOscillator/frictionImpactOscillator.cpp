/* Siconos-Examples , Copyright INRIA 2005-2011.
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
 * Contact: Maureen Chou maureenchou@.zju.edu.cn
*/
#include "SiconosKernel.hpp"
#include "frictionImpactOscillator.h"
using namespace std;
int main(int argc, char* argv[])
{
  boost::timer t;
  t.restart();
  try
  {
    // ================= Model definition =================

    // User-defined main parameters
    unsigned int nDof = 2;            // degrees of freedom
    double t0 = 0;                    // initial computation time
    double T = 3.0;                   // final computation time
    double h = 1e-5;               // time step
    double theta = 0.5;               // theta for MoreauJeanOSI integrator;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    SP::SiconosMatrix Mass;
    Mass.reset(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m_2;
    (*Mass)(1, 1) = J;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof));
    (*q0)(0) = y_0;
    (*q0)(1) = phi_0;

    SP::SiconosVector velocity0(new SiconosVector(nDof));
    (*velocity0)(0) = v_0;
    (*velocity0)(1) = omega_0;

    SP::LagrangianDS impactoscillator(new LagrangianDS(q0, v0, Mass));
    impactoscillator->setComputeFExtFunction("ImpactOsclillatorPlugin", "FExt");
    impactoscillator->setComputeFIntFunction("ImpactOsclillatorPlugin", "FInt");
    // --------------------
    // --- Interactions ---
    // --------------------
    string H = "ImpactOsclillatorPlugin:h";
	string jacH = "ImpactOsclillatorPlugin:jach";
	string dotH = "ImpactOsclillatorPlugin:doth";
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(eps_N, eps_T, mu, 2));
    SP::Relation relation(new LagrangianRheonomousR(H, jacH, dotH));
    SP::Interaction interaction(new Interaction(2, nslaw, relation));

    // -------------
    // --- Model ---
    // -------------

    SP::Model model(new Model(t0, T));
    model->nonSmoothDynamicalSystem()->insertDynamicalSystem(impactoscillator);
    model->nonSmoothDynamicalSystem()->link(interaction, impactoscillator);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretization --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
 
    SP::TimeStepping s(new TimeStepping(t));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator vOSI(new Moreau(impactoscillator, theta));
    s->insertIntegrator(vOSI);

    SP::OneStepNSProblem osnspb(new FrictionContact(2));
    s->insertNonSmoothProblem(osnspb);

    cout << "=== End of model loading === " << endl;

    // ================= Computation =================

    // --- Simulation initialization ---
    model->initialize(s);

    cout << "End of model initialisation" << endl;

    int k = 0;
    int N = floor((T - t0) / h);

    // --- Get the values to be plotted ---
    unsigned int outputSize = 5;
    SimpleMatrix dataPlot(N + 1, outputSize);
    dataPlot(k, 0) = t0;
    for (int i = 0; i < (int)nDof; i++)
    {
      dataPlot(k, 2 * i + 1) = (*impactoscillator->q())(i);
      dataPlot(k, 2 * i + 2) = (*impactoscillator->velocity())(i);
    }

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    while (k < N)
    {

      // get current time step
      k++;

      // solve ...
      s->computeOneStep();

      // get values
      dataPlot(k, 0) = s->nextTime();
      for (int i = 0; i < (int)nDof; i++)
      {
        dataPlot(k, 2 * i + 1) = (*dynamicalSystem->q())(i);
        dataPlot(k, 2 * i + 2) = (*dynamicalSystem->velocity())(i);
      }

      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
    }
    cout << "End of computation - Number of iterations done: " << k << endl;

    // --- Output files ---
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
    std::cout << "Comparison with a reference file" << std::endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("Woodpecker.ref", "ascii", dataPlotRef);
/*    double error = (dataPlot - dataPlotRef).normInf()/  dataPlotRef.normInf();
    std::cout << "Error = "<< error <<std::endl;
    if (error > 1e-12)
    {
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
      std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;
      return 1;
    }
*/
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "Exception caught in \'ImpactOsclillator\'" << endl;
    return 1;
  }
  cout << "Computation Time " << t.elapsed()  << endl;
  return 0;
}
