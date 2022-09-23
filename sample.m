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

% Download data u, dictionary matrix A, the legend of coefficient Xi, 
% and exact coefficient Xi_true.
load('aizawa_data.mat')
load('dictionary_matrix.mat')
load('dictionary_matrix_legend.mat')
load('aizawa_Xi_true.mat')
% add noise
eta = 1e-2;                      % variance
utilde = u+eta*randn(size(u));   % add noise
% numerical derivation
h = 0.05;                        % time step
b = derivative(utilde,h);
% find coefficients Xi
[n,d] = size(utilde); c = length(L);  % calculate Xi size
D = [5,6,7,8,9,10];                   % x_D may be dense
lambda0 = 0.09*ones(c,1);             % thresholding parameter
epsilon = 1;                          % auxiliary parameter
gamma = 0;                            % ridge parameter
Xi = zeros(c,d);
e = zeros(1,d);
% the iterative step
for i = 1:d
    x = RWPSTRidge(A,D,b(:,i),lambda0,epsilon,gamma);
    Xi(:,i) = x;  
    % relative error
    e(i) = norm(Xi(:,i)-Xi_true(:,i),2)/norm(Xi_true(:,i),2);
end

function udot = derivative(utilde,h)
% Calculate numerical derivative
[m,n] = size(utilde);                                        %calculate Xtilde size          
udot = zeros(m,n);                  
udot(1,:) = (utilde(2,:) - utilde(1,:))./h;                  
udot(2:m-1,:) = (utilde(3:m,:) - utilde(1:m-2,:))./(2*h);    %central difference method
udot(m,:) = (utilde(m,:) - utilde(m-1,:))./h;
end

