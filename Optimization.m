%% Steepest descent

clear
clc
close all

% Function and partial derivatives

% f =  @(x,y) (x-0.7).^2 + 3*(y-0.5).^2 + 0.5*x.^2.*y.^2;
% fx = @(x,y) 2*(x-0.7) + x.*y.^2;
% fy = @(x,y) 6*(y-0.5) + x.^2.*y;


% f = @(x,y) -cos(y.*x.^2) + (x-1).^2 .* (y-2).^2 + 1;
% fx = @(x,y) 2 * x .* y .* sin(y.*x.^2) + 2*(x-1) .* (y-2).^2;
% fy = @(x,y) x.^2 .* sin(y.*x.^2) + 2*(x-1).^2 .* (y-2);


% x = [1; -1] %Initial guess
x = [2; -2]

all_x = x;
iter_max = 50;
tol = 1e-4; % this is a fairly large tolerance... you may want to try with a smaller one
gamma = 1e-4; % parameter in Armijo line search

levels = [0.2:0.2:1, 1:6, 6:2:20];

[X,Y] = meshgrid(linspace(0,2,50),linspace(-2,2,50));
contour(X,Y,f(X,Y),levels,'b')
hold on
plot(all_x(1,:), all_x(2,:),'-or')
pause

%Initiliazation
iter = 0;
grad = [inf; inf]; % Gradient

while (norm(grad,2) >= tol) && iter < iter_max
    grad(1) = fx(x(1), x(2));
    grad(2) = fy(x(1), x(2));
    alpha = 1; %Step size
    x_new = x - alpha * grad;
    
    % Line search
    while f(x_new(1), x_new(2)) - f(x(1), x(2)) > - gamma * alpha * norm(grad,2)^2 && alpha > 10^-8
        alpha = alpha/2;
        x_new = x - alpha * grad;
    end
    
    clc
    alpha
    x = x_new
    f(x(1),x(2))
    all_x = [all_x x];
    plot(all_x(1,:), all_x(2,:),'-or')
    pause
    iter = iter + 1;
end

plot(all_x(1,:), all_x(2,:),'-or')
figure
plot(1:iter, f(all_x(1,2:end),all_x(2,2:end)), 'o')

%%

clear
clc
close all

% Rosenbrock's function, aka Steepest Descent's kryptonite. Google it!
f =  @(x,y) (1-x).^2 + 100*(y-x.^2).^2; % Min at (x,y) = (1,1)
fx = @(x,y) -2.*(1-x) - 400*x.*(y-x.^2);
fy = @(x,y) 200*(y-x.^2);

x = [0.5; 0.5]

all_x = x;
iter_max = 1000;
tol = 1e-4; % this is a fairly large tolerance... you may want to try with a smaller one
gamma = 1e-4; % parameter in Armijo line search

levels = [0:0.2:1, 1:0.5:5 6:11];

[X,Y] = meshgrid(linspace(0,1,50),linspace(0,1,50));
contour(X,Y,f(X,Y),levels,'b')
hold on
plot(all_x(1,:), all_x(2,:),'-or')
pause

%Initialization
iter = 0;
grad = [inf; inf]; % Gradient

while (norm(grad,2) >= tol) && iter < iter_max
    grad(1) = fx(x(1), x(2));
    grad(2) = fy(x(1), x(2));
    alpha = 1; %Step size
    x_new = x - alpha * grad;
    
    % Line search
    while f(x_new(1), x_new(2)) - f(x(1), x(2)) > - gamma * alpha * norm(grad,2)^2 && alpha > 10^-8
        alpha = alpha/2;
        x_new = x - alpha * grad;
    end
    
    clc
    alpha
    x = x_new
    f(x(1),x(2))
    all_x = [all_x x];
    plot(all_x(1,:), all_x(2,:),'-or')
    pause
    iter = iter + 1;
end

plot(all_x(1,:), all_x(2,:),'-or')
figure
plot(1:iter, f(all_x(1,2:end),all_x(2,2:end)), 'o')

%% Optimization in R
clear
clc

F = @(x) cos(5*x).*x.^2 - 1./((x-1).^2+0.1)
plot(linspace(0,2),F(linspace(0,2)));

% fmintx(F,0,2)
% fmintx(F,1,2)