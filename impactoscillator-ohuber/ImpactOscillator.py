
import Siconos.Numerics as SN
import Siconos.Kernel as SK
import numpy as np

import matplotlib.pyplot as plt

# system model parameters definition
m = 1.2
m_1 = 2.5
m_2 = 1.8
l = 1.5
Jt = 1.0 / 3 * m * l * l
ke = 2.0e4
Fl = 150.0
f_excite = 200
we = 2 * np.pi * f_excite
Fa = 50  # excitation force = Fa * sin(we * t)

# system contact parameters definition
eN = 0.2
eT = 0.0
mu = 0.1

# system initial state parameters definition
phi_0 = np.pi / 4.
omega_0 = 0.0
d0 =  l * np.sin(phi_0)
x1_0 = 0.0
v1_0 = 0.0
x2_0 = 0.0
v2_0 = 0.0

# system simulation parameters definition
nDof = 3  # degrees of freedom for robot arm
t0 = 0.         # initial computation time
T = 0.2       # final computation time
h = 1e-5       # time step : do not decrease, because of strong

# specifying initial conditions
q0 = np.zeros((nDof,))
# q0 --> x1, x2, phi
v0 = np.zeros((nDof,))
# v0 --> v1, v2, omega
q0[0] = x1_0
q0[1] = x2_0
q0[2] = phi_0

v0[0] = v1_0
v0[1] = v2_0
v0[2] = omega_0

    # ============================================
    # ============= Dynamical systems ============
    # ============================================
print("====> Model loading ...")

oscillator = SK.LagrangianDS(q0, v0, "ImpactOscillatorPlugin:mass")
oscillator.setComputeNNLFunction("ImpactOscillatorPlugin", "NNL")
oscillator.setComputeJacobianNNLqFunction("ImpactOscillatorPlugin", "jacNNLq")
oscillator.setComputeJacobianNNLqDotFunction("ImpactOscillatorPlugin", "jacNNLqDot")

oscillator.setComputeFIntFunction("ImpactOscillatorPlugin", "FInt")
oscillator.setComputeJacobianFIntqFunction("ImpactOscillatorPlugin", "jacFIntq")
oscillator.setComputeJacobianFIntqDotFunction("ImpactOscillatorPlugin", "jacFIntqDot")

oscillator.setComputeFExtFunction("ImpactOscillatorPlugin", "FExt")

# ============================================
# ================ Interactions ==============
# ============================================
nslaw = SK.NewtonImpactFrictionNSL(eN, eT, mu, 2)
relation = SK.LagrangianScleronomousR("ImpactOscillatorPlugin:Gap", "ImpactOscillatorPlugin:Gapjac")
inter =  SK.Interaction(2, nslaw, relation, 1)

# ============================================
# ================== Model ===================
# ============================================
impactOscillator = SK.Model(t0, T)
impactOscillator.nonSmoothDynamicalSystem().insertDynamicalSystem(oscillator)
impactOscillator.nonSmoothDynamicalSystem().link(inter, oscillator)

# ============================================
# ============== Simulation ==================
# ============================================
OSI = SK.MoreauJeanOSI(oscillator, 0.5, 0.0)
t = SK.TimeDiscretisation(t0, h)
osnspb = SK.FrictionContact(2, SN.SICONOS_FRICTION_2D_ENUM)

osnspb.numericsSolverOptions().dparam[0] = 1e-08
osnspb.numericsSolverOptions().iparam[0] = 100
osnspb.numericsSolverOptions().iparam[2] = 1; # random
s = SK.TimeStepping(t)
s.insertIntegrator(OSI)
s.insertNonSmoothProblem(osnspb, SK.SICONOS_OSNSP_TS_VELOCITY)
s.setNewtonTolerance(1e-10)
s.setNewtonMaxIteration(200)

topo = impactOscillator.nonSmoothDynamicalSystem().topology()
# ========= End of model definition=============

# ================================= Computation =================================

# ============================================
# ============== Computation =================
# ============================================
print("====> initialization ... It will soon be done!")
impactOscillator.initialize(s)
N = np.ceil((T - t0) / h) # Number of time steps

# ============================================
# --- Get the values to be plotted ---
# -> saved in a matrix dataPlot
# ============================================
outputSize = 11
dataPlot = np.empty((N + 1, outputSize))

q = oscillator.q()
v = oscillator.velocity()
y = inter.y(0);
lam = inter.lambda_(1)


# q :  x1 , x2 , phi
# v :  v1 , v2 , omega
dataPlot[0, 0] = impactOscillator.t0()
dataPlot[0, 1] = q[0]
dataPlot[0, 2] = q[1]
dataPlot[0, 3] = q[2]
dataPlot[0, 4] = v[0]
dataPlot[0, 5] = v[1]
dataPlot[0, 6] = v[2]

# --- Time loop ---

# ==== Simulation loop - Writing without explicit event handling =====
k = 1
while s.hasNextEvent() and k<= N:

    s.advanceToEvent()
    osnspb.setNumericsVerboseMode(0)
    # --- Get values to be plotted ---
    dataPlot[k, 0] = s.nextTime()
    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = q[1]
    dataPlot[k, 3] = q[2]
    dataPlot[k, 4] = v[0]
    dataPlot[k, 5] = v[1]
    dataPlot[k, 6] = v[2]
    dataPlot[k, 7] = y[0];
    dataPlot[k, 8] = y[1];
    dataPlot[k, 9] = lam[0];
    dataPlot[k, 10] = lam[1];



    s.processEvents()
    k += 1

print("End of computation - Number of iterations done: {:}".format(k - 1))

dataPlot.resize(k, outputSize)
np.savetxt("result-py.dat", dataPlot)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(dataPlot[:, 0], dataPlot[:, 1])
ax.set_xlabel('t')
ax.set_ylabel('q0')
fig.canvas.set_window_title('q0')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(dataPlot[:, 0], dataPlot[:, 2])
ax.set_xlabel('t')
ax.set_ylabel('q1')
fig.canvas.set_window_title('q1')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(dataPlot[:, 0], dataPlot[:, 3])
ax.set_xlabel('t')
ax.set_ylabel('q2')
fig.canvas.set_window_title('q2')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(dataPlot[:, 0], dataPlot[:, 4])
ax.set_xlabel('t')
ax.set_ylabel('v0')
fig.canvas.set_window_title('v0')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(dataPlot[:, 0], dataPlot[:, 5])
ax.set_xlabel('t')
ax.set_ylabel('v1')
fig.canvas.set_window_title('q1')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(dataPlot[:, 0], dataPlot[:, 6])
ax.set_xlabel('t')
ax.set_ylabel('v2')
fig.canvas.set_window_title('v2')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(dataPlot[:, 0], dataPlot[:, 10])
ax.set_xlabel('t')
ax.set_ylabel('lambda(1)')
fig.canvas.set_window_title('lambda(1)')



plt.show()
