clear
clc

f = @(x) (x.^2 + sin(x))./(1+x.^4);

plot(linspace(0,10), f(linspace(0,10)))

pause

quadgui(f, 0, 10, 1e-3)

quad(f, 0, 10, 1e-1)
quad(f, 0, 10, 1e-2)
quad(f, 0, 10, 1e-3)
quad(f, 0, 10, 1e-10)

%%

clear
clc

x  = linspace(-1,1,1000);
f = @(x) 1./(0.01+25*x.^2);

plot(x,f(x))

% the integral of f over [-1,1] equals 4*atan(50), which is about 6.2032
error = quad(f, -1, 1, 1e-6) - 4*atan(50) 

pause

quadgui(f, -1, 1, 1e-3)

%% Removable singularities

quadtx(@(x) sin(x)./x, 0, pi)

quadtx(@(x) sin(x)./x, realmin, pi)

quadtx(@sinc, 0, pi)

%% Integrals that depend on parameters

clear
clc

% Using a function handle
F = @(t,z,w) t.^(z-1).*(1-t).^(w-1);

z = 2;
w = 1.5;
tol = 1e-6;
quadtx(F,0,1,tol,z,w)

beta(z,w) % Matlab's built-in beta function