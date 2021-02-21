%% Backslash

A = [1 -1; 2 -1]

b = [2; 5]

x = A\b

%% Forward substitution (lower triangular systems)
clear
clc

n = 4;
L = tril(randi([1 5], n))  % a nxn lower triangular random matrix, whose elements are integers between 1 and 5
b = rand(n,1)
x = zeros(n,1);

pause
% solve Lx = b
for i = 1:n
    x(i) = b(i);
    for j = 1:i-1
        x(i) = x(i) - L(i,j)*x(j);
    end
    x(i) = x(i)/L(i,i);
end

x
pause

% Compare with backslash
y = L\b

%% Vectorize! You could save a lot of time.

clear
clc

n = 10000;
L = tril(rand(n));
b = rand(n,1);
x = zeros(n,1);

pause

tic % start timer

% solve Lx = b by doing forward substitution and a naive loop (as we just did)
x = zeros(n,1);
for i = 1:n
    x(i) = b(i);
    for j = 1:i-1
        x(i) = x(i) - L(i,j)*x(j);
    end
    x(i) = x(i)/L(i,i);
end

toc % stop timer

pause

tic
% solve Ly = b by doing a vectorized forward substitution 
y = zeros(n,1);
y(1) = b(1);
for i = 2:n
    y(i) = ( b(i) - L(i,1:i-1)*y(1:i-1) ) / L(i,i);
end
toc

%% Backward substitution (upper triangular systems)

clear
clc

n = 4;
U = triu(rand(n))
b = rand(n,1)
x = zeros(n,1);

x(n) = b(n)/U(n,n);
for i = n-1:-1:1
    x(i) = ( b(i) - U(i,i+1:n) * x(i+1:n) ) / U(i,i);
end

%% Permutation matrices
clc
clear

A = randi([0, 9],3)
P = [0 1 0; 1 0 0; 0 0 1]

pause

P*A % permutes rows
pause
A*P % permutes columns

% we could have done the same using a permutation vector
p = [2 1 3]
A(p,:) % permutes rows
A(:,p) % permutes columns

%% Product (and inverse) of two multiplier matrices
clear
clc

n = 5;
% $1 \le i, j < n$
i = 2; 
j = 4;

Mi = eye(n);
Mi(i+1:n, i) = -rand(n-i,1)

Mj = eye(n);
Mj(j+1:n, j) = -rand(n-j,1)

% inv_Mi = inv(Mi)
% Mi_times_Mj = Mi*Mj


%%  When solving lots of systems with the same (big) matrix A, do LU decomposition!

clear
clc

N = 100;
M = 100;

A = rand(N);
b = rand(N,M);
x = zeros(N,M);

tic
for i = 1:M
    x(:,i) = A\b(:,i);
end
lots_of_Gaussian = toc

tic
[L,U,P] = lu(A);

for i = 1:M
    y = L \ (P*b(:,i));
    x(:,i) = U \ y;
end
lu_once = toc

ratio = lots_of_Gaussian/lu_once



%%

clear
clc

A = [10 -7 0; -3 2 6; 5 -1 5]

[L,U,p] = lugui(A) % toolbox's interactive LU decomposition


%% Why is pivoting necessary? (computations from the example)
format long

% First step
2.099 + -7*0.3
3.901 + 7*0.3

% Second step
5 + 6 * 2.5e3
6.001 * 2.5e3

1.5002e4 + 2.5

% Backward substitution
% x3
(1.5004e4)/(1.5005e4)

% x2
6 * 0.99993
(6.001-5.9995) / (-1e-3)

% x1
-7*(-1.5)
(7-10.5)/10


%% A "bad" sysytem

A = [0.780 0.563; 0.913 0.659]
b = [0.217; 0.254]

%%
fig1 = ezplot(@(x,y) A(1,1)*x + A(1,2)*y - b(1), [-1.5, 1.5]);
hold on
fig2 = ezplot(@(x,y) A(2,1)*x + A(2,2)*y - b(2), [-1.5, 1.5]);
set(fig1, 'Color', 'r', 'linewidth', 2);
set(fig2, 'Color', 'b', 'linewidth', 2);


A(1,1)/A(2,1)
A(1,2)/A(2,2)
b(1)/b(2)


%% Condition number

cond(A,1)


db = 10^-5*[-1; 1]


x1 = A\b
x2 = A\(b+db)

norm(db)/norm(b)
norm(x2-x1)/norm(x1)

%% Matrix from Homework #2, exercise 7
n = 40;
A = diag(2 * ones(n, 1), 0) - diag(ones(n-1, 1), -1) - diag(ones(n-1, 1), 1);

nnz(A)
density = nnz(A)/prod(size(A))
sparsity = 1 - density


%%
B = sparse(A);

% use whos to compare how much space A and B take up, especially for n
% large

%%
clear
clc

i = [1 1 2 4];
j = [1 3 1 2];
x = [5 5 5 5];

A = sparse(i,j,x,4,4)
B = full(A)

%%
% Bandwidth
[i,j] = find(A);
bandwith = max(abs(i-j))

%%
clear
clc

n = 4;
a = - ones(n-1,1);
b = 2* ones(n,1);
c = a;
d = [1:n];

tridisolve(a,b,c,d)