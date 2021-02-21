%% Solving the ODE y'(t) = a*y(t), y(0) = 1. The exact solution is y(t) = e^(a*t).

clear 
clc
close

a = -10;
T = 10;
n = 20;
h = T/n;
y = zeros(1,n+1);
y(1) = 1;

%%
% Forward Euler
for j = 1:n
    y(j+1) = y(j) * (1 + a*h);
end

plot(linspace(0,T,n+1),y, '-b', 'LineWidth', 2)
hold on
plot(linspace(0,T,n+1),exp(a*linspace(0,T,n+1)), '-r', 'LineWidth', 2)


%%
close

% Backward Euler
for j = 1:n
    y(j+1) = y(j) / (1 - a*h);
end

plot(linspace(0,T,n+1),y, '-b', 'LineWidth', 2)
hold on
plot(linspace(0,T,n+1),exp(a*linspace(0,T,n+1)), '-r', 'LineWidth', 2)

%% Matlab solvers

a = -10;
F = @(t,y) a*y;
tspan = [0, 1];
y0 = 1;

% opts = odeset('reltol',1e-3,'abstol', 1e-5);
% [tout, yout] =  ode23(F,tspan,y0);
[tout, yout] =  ode23(F,tspan,y0);

plot(tout, yout, 'o-')
hold on
plot(linspace(0,1), y0*exp(a*linspace(0,1)));


%%

F = @(t,y) (1-t)./(1+t.*y^2);
tspan = [0, 10];
y0 = 1;
% ode23tx(F,tspan,y0)

ode23tx(F,tspan,y0);
% plot(tout, yout, 'o-')


%% Harmonic oscillator (mass and spring)

k = 1;
mu = -0.1; % friction
F = @(t,y) [0 1; -k mu]*y;
tspan = [0, 20];
y0 = [0; 1];

[tout, yout] = ode23(F,tspan,y0,1e-6);
% plot(yout(:,1),yout(:,2),'-o')
plot(tout,yout(:,1),'b', tout,yout(:,2),'r')

%% Stiffness (lighting a match, see section 7.9)

clear
clc

F = @(t,y) y^2 - y^3;
delta = 0.00001;
y0 = delta;
tspan = [0 2/delta];
opts = odeset('RelTol',1.e-4);
ode23s(F, tspan, y0, opts);

% [tout, yout] = ode23s(F,tspan,y0,1e-6);
% plot(tout,yout,'-o')

%% Events (Falling body)

clear 
clc

F = @(t,y) [y(2); -1+y(2)^2];
tspan = [0, 5];
y0 = [1; 0];
[tout, yout] = ode23(F,tspan,y0);
plot(tout, yout(:,1), 'o-', tout, 0*tout, '-k')

%% 

clear
clc

F = @(t,y) [y(2); y(2)^2-1];
tspan = [0, 20];
y0 = [1; 0];
opts = odeset('events',@g_fall);
[tout, yout, t_final] = ode23(F,tspan,y0,opts);

plot(tout, yout(:,1), '-o')

t_final
%% Pendulum

clear 
clc

g = 9.81;
L = 1;
alpha = 0.;
F = @(t,y) [y(2); -g/L*sin(y(1))-alpha*y(2)];
tspan = [0, 10];
y0 = [pi/4; 0];
tol = 1e-3;

[tout, yout] = ode113(F,tspan,y0,tol);
% length(tout)
plot(tout, yout(:,1),'o-')
hold on
plot(tout, yout(:,2),'or-')
% plot(yout(:,1), yout(:,2),'o-')

%%

clear 
clc

g = 9.81;
L = 2;
alpha = 0;
F = @(t,y) [y(2); -g/L*sin(y(1))-alpha*y(2)];
tspan = [0, 10];
y0 = [pi/4; 0];

opts = odeset('events', @(t,y)g_pendulum(t,y,y0), 'reltol', 1e-4);
[tout, yout, t_final, y_final] = ode23(@(t,y)F(t,y),tspan,y0,opts);
t_final
plot(tout, yout(:,1),'o-',tout, yout(:,2),'o-r')