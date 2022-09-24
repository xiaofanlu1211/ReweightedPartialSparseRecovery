% =========================================================================
% Goal: Find the value of the objective function F_{wp}(x)
%       and the change of the threshold parameters \lambda
% 
% Authors: Xiaofan Lu, Huimei Ma, and Linan Zhang
% Date: Sept, 2022
% =========================================================================
% Construct matrix A and vector b
m = 10; n = 10;                       % size of matrix A       
A = randi(n,m,n);                     % construct an m-by-n matrix of pseudorandom integers
x_true = zeros(n,1);                  % construct exact vector x
x_true(1:2:5) = 1;
eta = 0.1*randn(m,1);                 % construct noise eta~N(0,0.1^2)
b = A*x_true + eta;                   % construct vector b:= Ax+eta
D = linspace(1,2,2);                  % dense part/prior known support set of x
lambda0 = 0.1*ones(n,1);              % initial threshold parameter 
epsilon = 1;                          % reweighted auxiliary parameter
gamma = 0;                            % ridge parameter
% Find sparse approximation solution of x
% Find the objective function F_{wp}(x)
[x,lambda,Xout,Xin,kout,kin,F] = RWPSTRidge(A,D,b,lambda0,epsilon,gamma);
% compute relative error
 e = zeros(kin+1+kout,1);
 for ii =1:kin+1+kout
     e(ii) = norm(Xin(:,ii) - x_true,2)/ norm(x_true,2);
 end
