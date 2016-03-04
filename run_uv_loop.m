% function run_uv_loop

% This script solves U_t = F(U,t)  using a pseudospectral method.
%==========================================================================
ord = 3; % n = 2^ii  space grid points, increase for more accuracy
n_loop = 3; %how many loops do you want to run over?



for ii = ord:1:ord +n_loop %loop for n_loops at resolution ii
%% Set options
%%
plotuv = 'yes';      % 'yes' to plot u and v
plotha = 'yes';     % 'yes' to plot H1 and A norms
holdon = 'yes';     % 'yes' to keep old H1 and A-norm plots
lineha = 'r:';     % linestyle for norm plots

%% Set space parameters
%%
display(ii)

x0 = 0.;             % left endpoint
x1 = 2.*pi;          % right endpoint
xlength = x1-x0;     % length of interval solving over 
n = 2^ii;            % number of points
dx = xlength/n;      % distance between points on grid
x = x0+dx*(0:(n-1)); % create spatial grid
lx = length(x);      % 
k = make_k(lx);

%% Set time parameters
%%
t0 = 0.;             % initial time
tf = 1.;             % final time
dt=2^(-8);
ts = 2^3;
t= t0:dt:tf;
lt = length(t);





%% Set coefficients
%%
alpha = 1; %coefficient of nonlinear term
beta = 2; %coefficient of nonlinear term
nu = 0; % coefficient of dispersive term

%% Define initial condition  uv0 = uv(x,0)
%%
amp1 = 1.;  amp2 = 2.;
u0 = amp1.*exp(1i.*x) + amp2.*exp(2i.*(x+2*pi^2));
v0 = amp1.*exp(-1i.*x) + amp2.*exp(-2i.*(x+2*pi^2));
u0 = P(u0,lx);
u0 = truncate(u0,lx);
v0 = 0*Q(v0,lx);
v0 = truncate(v0,lx);
uv0 = [u0 v0].';

switch plotuv
    case 'yes'
        uplot = zeros(lt,lx);   %stores values of u
        vplot = zeros(lt,lx);   %stores values of v
        uplot(1,:) = u0;    %input initial value of u
        vplot(1,:) = v0;    %input initial value of v
end

switch plotha
    case 'yes'
        h1 = t;
        anorm = t; 
        upv=u0+v0;
        h1(1) = sum(abs(deriv(upv,k)).^2)/lx;   %store initial h1 norm
        anorm(1) =  sum(abs(fft(upv)))/lx;  %store initial a norm     
end

%% Solve the equation
%%
tj=t0;
ntcount = 0;
tic
for jj=2:2:lt-2 %cycle through time span
    tspan = [t(jj), t(jj+1), t(jj+2)]; %generate time span to solve on
    [~, uv] = ode45(@f_uv,tspan,uv0,[],k,lx,alpha,beta,nu); %solve a N system of ODE's
    for kk = 2:3 %plot solutoin at tj+ dtplot, +2 dtplot
        u0(:) = uv(kk,1:lx); %store solution
        v0(:) = uv(kk,lx+1:2*lx);
        
        switch plotuv
            case 'yes'
                uplot(jj+kk-2,:) = u0;
                vplot(jj+kk-2,:) = v0;
        end
        
        switch plotha
            case 'yes'        
                upv = u0 + v0;
                h1(jj+kk-2) = sum(abs(deriv(upv,k)).^2)/lx;
                anorm(jj+kk-2) =  sum(abs(fft(upv)))/lx;
        end  
        
        uv0 = [u0 v0].';
        ntcount = ntcount+1;
        if ntcount >= ntskip
            display(t(jj))
            ntcount=0;
        end
    end   
end
toc

%% Plot the solution
%
%
switch plotuv
    case 'yes'
        moduplot = abs(uplot + vplot);
 
        h = figure(71);  clf;

    tic
        waterfall(x,t,moduplot)
    toc   
        view(10,70);
        axis tight
        
        xlabel('x'), ylabel('\tau'), zlabel('|A|')
        grid off
        drawnow


        
        
% tic
%   name = [ 'figure' int2str(ii)  ];
% 	print(h,'-depsc','-r300',['fancy' name '.eps' ])
% toc
end

switch plotha
    case 'yes'
        figure(72);              % H^1-norm int[|D\psi|] vs. t
        switch holdon
            case 'yes'
                hold on ;
            otherwise
                clf;
        end
        emax=max(0,1.2*max(h1)); emin = min(0,1.2*min(h1));
        plot(t,h1, lineha);
        axis ([t0 tf emin emax]);
        xlabel('\tau'); ylabel('H^1-norm');

        figure(73);               % A-norm int[|\hat{\psi}|] vs. t 
        switch holdon
            case 'yes'
                hold on ;
            otherwise
                clf;
        end
        emax=max(0,1.2*max(anorm)); emin = min(0,1.2*min(anorm));
        plot(t,anorm,lineha);
        axis ([t0 tf emin emax]);
        xlabel('\tau'); ylabel('A-norm');
end
% end


  




end




