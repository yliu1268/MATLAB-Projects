%%
n = 20;
x = linspace(0.01,1,n+1);
y = exp(x.*sin(x.^2)) - 0.5*cos(1./x);
plot(x,y,'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'k')


%%
subplot(2,2,1)
plot(x,y,'ok', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
hold on
vector = linspace(min(x),max(x),200);
pol = polyinterp(x,y,vector); % interpolating polynomial
plot(vector,pol,'-g', 'LineWidth', 2)

pause

subplot(2,2,2)
plot(x,y,'ok', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
hold on
lin = piecelin(x,y,vector); % piecewise linear
plot(vector,lin,'-k', 'LineWidth', 2)

pause

subplot(2,2,3)
plot(x,y,'ok', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
hold on
spline_fun = splinetx(x,y,vector); % piecewise cubic spline
plot(vector,spline_fun,'-r', 'LineWidth', 2)

pause

subplot(2,2,4)
plot(x,y,'ok', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
hold on
shape_preserving = pchiptx(x,y,vector); % shape-preserving piecewise cubic
plot(vector,shape_preserving,'-b', 'LineWidth', 2)

%% Vandermonde matrices are ill-conditioned

clear
clc
n = 10;
x = linspace(0,1,n)';
cond(vander(x))

%% Interpolating polynomial and Runge's function

clear
clc

n = 20;
x = linspace(-1,1,n)';
% y = @(z) (1-z).*sin(pi*z)+z.^2.*cos(z); 
y = @(z) 1./(1+25*z.^2);

vector = linspace(-1,1,1000);
pol = polyinterp(x,y(x),vector);

% Plot
plot(x,y(x),'ok', 'MarkerSize', 8,'MarkerFaceColor', 'k')
hold on
plot(vector, y(vector),'-k', 'LineWidth', 2);
hold on 
pause
plot(vector,pol,'-r', 'LineWidth', 2)

%% Piecewise linear interpolation

clear
clc

n = 21;
x = linspace(-1,1,n)';
y = @(z) 1./(1+25*z.^2);

vector = linspace(-1,1,1000)';
pol = polyinterp(x,y(x),vector); % interpolating polynomial
lin = piecelin(x,y(x),vector); % piecewise linear interpolation

% Plot
plot(x,y(x),'ok', 'MarkerSize', 8,'MarkerFaceColor', 'k')
hold on
plot(vector, y(vector),'-k', 'LineWidth', 2);
hold on 
pause
plot(vector,lin,'-b', 'LineWidth', 2)
% pause
% plot(vector,pol,'-r', 'LineWidth', 2)

% Errors?
% error_pol = max(abs(y(vector)-pol))
% error_lin = max(abs(y(vector)-lin))

%% Piecewise linear interpolation

clear
clc

n = 6;
x = linspace(0,2,n)';
% y = @(z) exp(2*z);

y = @(z) sin(z);

vector = linspace(0,2,1000)';
lin = piecelin(x,y(x),vector); % piecewise linear interpolation

% Plot
% plot(x,y(x),'ok', 'MarkerSize', 8,'MarkerFaceColor', 'k')
% hold on
% plot(vector, y(vector),'-k', 'LineWidth', 2);
% hold on 
% plot(vector,lin,'-b', 'LineWidth', 2)

error_lin = max(abs(y(vector)-lin))

%%

x = [0 1];
y = [0 4];

v = piecelin(x,y, [0.1, 0.7])


%% Errors: piecewise cubic Hermite vs. piecewise linear

% clear
% clc

n = 12;
x = linspace(0,1,n);
f = @(z) exp(2*z);
f_prime = @(z) 2*exp(2*z);

% Try with any function and its derivative
% f = @(z) (1-z).^2.*sin(10*z);
% f_prime = @(z) -2*(1-z).*sin(10*z)+10*(1-z).^2.*cos(10*z);

u = linspace(0,1,1000);
v = linspace(0,1,1000);
index = 1;

for k = 1:n-1
    % Find Hermite interpolant at every subinterval. Here, we are writing the
    % piecewise cubic polynomial as 
    % p(x) = y(x(k)) + y_prime(x(k))*(x-x(k)) + a*(x-x(k))^2 + b*(x-x(k))^3
    % and finding a and b (using the information we have on x(k+1)). 
    % Note that this is not the formula we used in class: see line 33 in
    % either splinetx or pchiptx.
    
    delta_x = x(k+1)-x(k);
    M = [delta_x^2 delta_x^3;
         2*delta_x 3*delta_x^2];
    b = [f(x(k+1))-f(x(k)) - f_prime(x(k))*delta_x; f_prime(x(k+1))-f_prime(x(k))];
    sol = M\b;
    hermite_pol = @(z) f(x(k)) + f_prime(x(k))*(z-x(k)) + sol(1)*(z-x(k))^2 + sol(2)*(z-x(k))^3;
    % Evaluate polynomial at u
    while u(index) <= x(k+1) && index < 1000
        v(index) = hermite_pol(u(index));
        index = index + 1;
    end
    v(end) = f(x(end));
end

piecewise_linear = piecelin(x,f(x),u);


% Plot
% plot(x,f(x),'ok', 'MarkerSize', 8)
% hold on
% plot(u, f(u),'k', 'LineWidth', 2);
% hold on 
% pause
% plot(u,v,'r', 'LineWidth', 2)
% hold on
% plot(u,piecewise_linear,'b', 'LineWidth', 2)

err_hermite = max(abs(f(u)-v)) % Piecewise cubic error 
err_linear = max(abs(f(u)-piecewise_linear)) % Piecewise linear error

%%
n = 20;
x = linspace(0,1,n);
y = cos(exp(x.^2)) + 0.5*x.*rand(1,n).^2;
plot(x,y,'ok', 'MarkerSize', 8, 'MarkerFaceColor', 'k')

pause

vector = linspace(min(x),max(x),200);

subplot(1,2,1)
plot(x,y,'ok', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
hold on
spline_fun = splinetx(x,y,vector); % piecewise cubic spline
plot(vector,spline_fun,'-r', 'LineWidth', 2)

pause

subplot(1,2,2)
plot(x,y,'ok', 'MarkerSize', 6, 'MarkerFaceColor', 'k')
hold on
shape_preserving = pchiptx(x,y,vector); % shape-preserving piecewise cubic
plot(vector,shape_preserving,'-b', 'LineWidth', 2)