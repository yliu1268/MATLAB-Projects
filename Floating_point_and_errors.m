%% What's wrong with MATLAB?

x = 0.5/0.1 - 5

y = 0.3/0.1 - 3


%% Solving an overdetermined system?

A = [17 6; 1.7 0.6];
b = [23; 2.3];
x = A\b

%% What about this?

small_number = 1e-16;
1 + small_number - 1
1 - 1 + small_number

%%

x = 0.988:.0001:1.012;
% y = x.^7-7*x.^6+21*x.^5-35*x.^4+35*x.^3-21*x.^2+7*x-1; % note that y = (x-1).^7
% plot(x,y)
z = (x-1).^7;
plot(x,z)

%% Get the book's toolbox!
% t: # bits in the mantissa
% [emin,emax]: range for the exponent

floatgui

%% Some special numbers

eps
realmin
realmax
2*realmax
Inf-Inf
0/0
eps*realmin
2.47e-324
2.48e-324

%% Roundoff error
% The colon command avoids propagation of roundoff errors
clc

x = 0;
t = 0.1;

for i = 1:10
    x = x + t;
end

1-x % x is not equal to 1!

y = 0:0.1:1;

1-y(end)

%% Approximation error
% Want to approximate $\int_0^1 e^{\sin x} dx$.

clear
clc

step = 0.1;
subintervals = 0:step:1;
sum = 0;

for i = 1:size(subintervals,2)-1
    sum = sum + exp( sin(subintervals(i))) * step;
end

sum % this is our Riemann sum that approximates $\int_0^1 e^{\sin x} dx$.

%% Error propagation

x = sqrt(3);
y = x.^2 - 1; % y should equal 2
eps_y = (y-2)/2 % relative error in y; its absolute value should be less than 5*eps 

%% Numerical instability: we should have I_n -> 0 as n -> infinity

% We want to compute I_n = quad(@(x) x.^n./(x+5),0,1) using a recursive
% formula

p = zeros(1,10);
p(1) = log(6/5); % p(1) = I_0
for n = 1:length(p)
    p(n+1) = 1/(n+1) - 5*p(n);
end
% p

pause

p(10)
quad(@(x) x.^9./(x+5),0,1) % this is what p(10) should be!

%% Same thing, although now we move backwards in n

q = zeros(1,30);
q(end) = 0; % just a guess... anyway, I_30 should be very small!

for n = length(q)-1 :-1 :1 % note that n is moving backwards
    q(n) = (1/(n+1) - q(n+1))/5;
end
% q

q(9) % should equal quad(@(x) x.^9./(x+5),0,1)