# SPP_Simulation
Numerical Simulation of Surface Plasmon Polaritons using spectral and time stepping methods.


Main Scripts

run_uv_loop.m - main file, run this to find solution to PDE with initial defined within

f_uv.m - spectral method redues PDE to ODE of form u_t = f(u,t)

ode4_step.m - implementation of runge-kutta 4th order to solve ODE

