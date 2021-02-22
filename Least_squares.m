%% Least squares: motivation

clear 
clc
close 

err = 5*10^-2; % some measurement error
n = 10;

t = linspace(0,1,n)'; 
y = (t-t.^2) .* (1 + err * (2*rand(1,n)' - 1)); % parabola with some noise
% y = rand(1,n)';
plot(t, y, '*r')

pause 

u = linspace(0,1); % vector to make the plot look nice
a = polyfit(t,y,2);
hold on
plot(u, a(1)*u.^2+a(2)*u + a(3), 'b')


% Other way to compute this parabola: using backslash
% Basis: phi_1(t) = t^2, phi_2(t) = t, phi_3 (t) = 1

X = [t.^2, t, ones(n,1)];
beta = X\y;

%% Weighted least squares: multiply the design matrix and the observations by the weight vector.

clear 
clc
close 


err = 5*10^-2;
n = 10;

t = linspace(0,1,n)'; 
y = (t-t.^2) .* (1 + err * (2*rand(1,n)' - 1));
plot(t, y, '*r')
hold on

% same example as above, but now suppose we only trust observations 
% 1, 3, 4 and 6. In this example, we are going to discard the rest of
% the observations (setting w_i = 0 for i different from {1,3,4,6}); we
% could also give the other observations a different weight.

w = zeros(1,n);
w(1) = 1;
w(3) = 1;
w(4) = 1;
w(6) = 2;

pause 

X = [t.^2, t, ones(n,1)];
beta = (diag(w)*X)\(diag(w)*y);

u = linspace(0,1); % vector to make the plot look nice
plot(u, beta(1)*u.^2+beta(2)*u + beta(3), 'b')


%% Householder transformations: compare this with qrsteps and see end of section 5.5 in Moler's book.

clear
clc

s = (0:0.2:1)' % a 6x1 vector
X = [s.^2, s, ones(size(s))] % 6x3 matrix: we are using quadratic polynomial as a model

sigma = norm(X(:,1)) % X(1,1) = 0, we don't need to worry about the sign of sigma

u = X(:,1) + [sigma; 0; 0; 0; 0; 0]

rho = 2/(u'*u)

H1 = eye(6) - rho*u*u' % first Householder transformation

X1 = H1*X % first step

pause

clc

X2_tilde = X1(2:end,2)

sigma = norm(X2_tilde(:,1))*sign(X2_tilde(1,1))

u = X2_tilde(:,1) + [sigma; 0; 0; 0; 0]

rho = 2/(u'*u);

H2_tilde = eye(5) - rho*u*u'

H2 = [1 0 0 0 0 0; zeros(5,1) H2_tilde] % second Householder transformation

X2 = H2*X1 % second step

pause

clc

X3_tilde = X2(3:end,3)

sigma = norm(X3_tilde(:,1))*sign(X3_tilde(1,1))

u = X3_tilde(:,1) + [sigma; 0; 0; 0]

rho = 2/(u'*u);

H3_tilde = eye(4) - rho*u*u'

H3 = [1 0 0 0 0 0; 0 1 0 0 0 0; zeros(4,2) H3_tilde] % third Householder transformation

X3 = H3*X2 % third step


%% Resulting matrices

clc

R = X3
Q = H1*H2*H3 % in practice, we won't need to compute this matrix

% Check that Q is orthogonal
% Q'*Q
% Q*Q'

%% Use this in an application: example from Section 5.5

y = [150.697;
     179.323;
     203.212;
     226.505;
     249.633;
     281.422]
 
 z = H3*H2*H1*y % 6x1 vector. 
 % The first 3 elements in z are related to the least square solution beta,
 % while the last 3 elements are the residual.
 
 R_tilde = R(1:3,1:3); % get rid of the elements below the diagonal (recall they are all 0)
 beta = R_tilde\z(1:3) % R_tilde is upper triangular: we can use backward substitution
 
 residual = z(4:6)
 norm(residual,2)
 
%% This is what Matlab's backslash does in this case!
bslash_sol = X\y % equals beta

bslash_residual = y - X*bslash_sol % not equal to our vector residual, but...?

%% Pseudoinverse of a matrix

clear all
clc

X = [1 1; 1 1]
Z = pinv(X) % Moore-Penrose pseudoinverse

norm(X*Z-eye(2), 'fro') % compute Frobenius norm
norm(Z, 'fro')

pause

% another minimizing matrix
U = [1/2 1/2; 0 0]
norm(X*U-eye(2), 'fro') % this equals to norm(X*Z-eye(2), 'fro')
norm(U, 'fro') % this is larger than norm(Z, 'fro')

%% Least squares with a rank-deficient design matrix

clear all
clc

M = randi(5,6,2);
X = [M M(:,1)+ M(:,2)]

rank(X) % X is a rank-two matix

y = rand(6,1);

pause 

beta_backslash = X\y 
beta_pinv = pinv(X)*y 

pause
norm(beta_backslash)
norm(beta_pinv)