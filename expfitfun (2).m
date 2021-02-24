   function res = expfitfun(lambda,t,y)

   % Computes the solution to the linear least square that results when the
   % lambdas are fixed.
   % res equals the norm of the residual 
   
   m = length(t);
   n = length(lambda);
   X = zeros(m,n);
   
   % Design matrix (using the lambdas we have)
   for j = 1:n
      X(:,j) = exp(-lambda(j)*t);
   end
   
   % Solve a linear LS problem and compute the residual.
   beta = X\y;
   z = X*beta;
   res = norm(z-y);
   
   end