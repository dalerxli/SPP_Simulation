# SPP_Simulation
Numerical Simulation of Surface Plasmon Polaritons using spectral and time stepping methods.


Main Scripts

run_uv_loop.m - Main file, run this to find solution to PDE with initial defined within

f_uv.m - Spectral method redues PDE to ODE of form u_t = f(u,t), f_uv.m gives right hand side

ode4_step.m - Implementation of 4th order Runge-Kutta to solve ODE


Remaining scripts are implementations of spatial derivatives and operators used within the main scripts.
