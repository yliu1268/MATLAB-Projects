% Strongly recommended: read section 5.8
% Here we are fitting to a model y(t) = beta_1 * exp(-lambda_1*t) + beta_2 * exp(-lambda_2*t)

clear
clc

load data.txt

t = data(:,1);
y = data(:,2);

lambda0 = [-1 -2]'; % initial guess
m = length(t); % number of data points
n = length(lambda0);

fun = @(lambda) expfitfun(lambda,t,y);

lambda = fminsearch(fun,lambda0) % find the optimal lambdas
   
 X = zeros(m,n); 
 for j = 1:n
     X(:,j) = exp(-lambda(j)*t);
 end
 
beta = X\y % solve a linear LS problem to find the betas
z = X*beta;
res = norm(z-y);

plot(t,y,'o');
hold on
vec = linspace(min(t),max(t));
plot(vec, beta(1)*exp(-lambda(1)*vec)+beta(2)*exp(-lambda(2)*vec));