function [x,lambda,Xout,Xin,kout,kin,F] = RWPSTRidge(A,D,b,lambda0,epsilon,gamma)
% =========================================================================
% Goal: to solve the linear system Ax=b using the proposed
%       algorithm when partial support set D of x known a- 
%       priori.
%
% Inputs:
%   A = matrix of size m*n
%   D = expect dense part
%   b = vector of size m*1
%   lambda0 = initial thresholding parameter
%   epsilon = auxiliary parameter of update lambda
%   gamma = ridge parameter
%
% Outputs:
%   x = the approximate solution of linear system Ax=b
%   lambda = the lambda at each iteration
%   Xout = the solution at each out loop iteration
%   Xin = the solution at each inner loop iteration
%   kout = iteration times of outer loop
%   kin = total iteration times of inner loop
%   F = the value of objective function F_{wp}(x)
% Authors: Xiaofan Lu, Huimei Ma, and Linan Zhang
% Date: Thursday, Sept 1, 2022
% =========================================================================

% parameters
n = length(lambda0); err = 1;
MaxIt = n;                     % Maximum number of iterations                    
tol = 1e-6;                    % tolerance
kout = 0;                      % number of outer iterations
kin = 0;                       % total number of inner iterations
Xout = zeros(n,MaxIt);         % record x^{k_out}
lambda = zeros(n,MaxIt);       % record \lambda^{k_out}
lambda(:,1) = lambda0;
% the iterative step
while kout <= MaxIt && err>tol
    % solve the linear system Ax=b
    [E,X,x,K] = WPSTRidge(A,D,b,lambda(:,kout+1),gamma);
    % update parameters, matrix and compute relative error
    Xout(:,kout+1) = x;
    if kout == 0
       Xin = X;
       F = E;
    else
       Xin = [Xin X];
       F = [F; E];
       err = norm(Xout(:,kout+1)-Xout(:,kout))/norm(Xout(:,kout+1));
    end    
    
    kout = kout + 1;
    kin = kin + K;
    lambda(:,kout+1) = lambda0.*(1./(abs(x)+epsilon));
end
Xout(:,kout+1:end) = [];
lambda(:,kout+1:end) = [];
kout = kout - 1;
end