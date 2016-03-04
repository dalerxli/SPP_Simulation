function y = ode4_step(odefun,ti,hi,yi,k,disp,damp,x,lx,alpha,beta)

%ODE4_STEP implements a 4th order Runge-Kutta method
%ti = initial time, hi = time step, yi = initial data


neq = length(yi);

F = zeros(neq,4);

F(:,1) = feval(odefun,ti,yi,k,disp,damp,x,lx,alpha,beta);

F(:,2) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,1),k,disp,damp,x,lx,alpha,beta);

F(:,3) = feval(odefun,ti+0.5*hi,yi+0.5*hi*F(:,2),k,disp,damp,x,lx,alpha,beta);  

F(:,4) = feval(odefun,ti+hi,yi+hi*F(:,3),k,disp,damp,x,lx,alpha,beta);

y = yi + (hi/6)*(F(:,1) + 2*F(:,2) + 2*F(:,3) + F(:,4));



