% =========================================================================
% Target equation: the Aizawa equations
%      udot_{1} = (u_{3}-0.7)*u_{1} - 3.5*u_{2};
%      udot_{2} = 3.5*u_{1} + (u_{3}-0.7)*u_{2};
%      udot_{3} = 0.6+0.95*u_{3} - u_{3}^3/3 -...
%                 (u_{1}^2+u_{2}^2)*(1+0.25*u_{3}) + 0.1*u_{3}*u_{1}^3
%
% Auxiliary functions:
%   - library.m (The Candidate Functions)
%   - derivative.m (Numerical derivation) 
%   - RWPSTRidge.m (Reweighted Partial STRidge Algorithm)
%   - aizawa.m (ODE to Test)
%
% Authors: Xiaofan Lu, Huimei Ma, and Linan Zhang
% Date: Aug, 2022
% =========================================================================

% =========================================================================
% Steps 1-4: Construct dictionary matrix A and velocity vector b.
% =========================================================================
% Step 1: make data
h = 0.05;                        % time step
u = aizawa(h);                   % the true Aizawa system
X = u';                                 
% Step 2: add noise
eta = 1e-2;                      % variance
Xtilde = X+eta*randn(size(X));   % add noise
% Step 3: numerical derivation
b = derivative(Xtilde,h);
% Step 4: construct dictionary matrix
% library of polynomials of degree at most p
% library of known polynomials of degree s
p = 5; s = 2;
[D,A,L] = library(Xtilde,p,s);
% Step 5: find coefficients Xi
[n,d] = size(Xtilde);           % calculate Xi size
c = length(L);
lambda0 = 0.09*ones(c,1);       % thresholding parameter
epsilon = 1;                    % auxiliary parameter
gamma = 0;                      % ridge parameter
Xi = zeros(c,d);
e = zeros(1,d);
% exact coefficient Xi
load ('Xi_true.mat')
% the iterative step
for i = 1:d
    x = RWPSTRidge(A,D,b(:,i),lambda0,epsilon,gamma);
    Xi(:,i) = x;
    % relative error
    e(i) = norm(Xi(:,i)-Xi_true(:,i),2)/norm(Xi_true(:,i),2);
end


