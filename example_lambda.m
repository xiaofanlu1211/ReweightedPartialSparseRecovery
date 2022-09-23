% =========================================================================
% Goal: Test whether RWPSTRidge algorithm is less dependent choice of 
%       lambda than WPSTRidge algorithm.

% Authors: Xiaofan Lu, Huimei Ma, and Linan Zhang
% Date: Sept, 2022
% =========================================================================
% Construct matrix A and vector b
m = 8; n = 8;                % size of matrix A
A = randn(m,n);              % Generate a m-by-n matrix of normally distributed random numbers
A = tril(A);                 % Lower triangular matrix   
A = normc(A);                % normalize the columns of A
x = zeros(n,1);              % ''true'' vector x
x(1:3) = 1;
eta = randn(m,1)*0.01;       % random noise eta~N(0,0.01^2)
b = A*x + eta;               % construct vector b:= Ax+eta
D = linspace(1,4,4);         % dense part/prior known support set of x
lambda0 = 0.1*ones(n,1);     % initial threshold parameter
epsilon = 1;                 % reweighted auxiliary parameter
gamma = 0;                   % ridge parameter

% Apply the RWPSTRidge algorithm to solve Ax=b
[x,lambda,Xout,Xin] = RWPSTRidge(A,D,b,lambda0,epsilon,gamma);