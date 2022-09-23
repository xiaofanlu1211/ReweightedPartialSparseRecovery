function u = aizawa(h)
% ============================================================
% Description:
%   The function outputs data u
%
% Inputs:
%   h = time step
%
% Outputs:
%   u = u(t), t = t_1,t_2,...,t_n
% Authors: Linan Zhang
% ============================================================

% Step 1: Set parameters.
T = 250;
t = 0:h:T;

% Step 2: Initialize the solution.
n = length(t);
u = zeros(3,n);
u(:,1) = [1; 1; 0.5];

% Step 3: Define u', denoted by f.
f = @(u) [(u(3)-0.7)*u(1) - 3.5*u(2);
    3.5*u(1) + (u(3)-0.7)*u(2);
    0.6+0.95*u(3) - u(3)^3/3 - (u(1)^2+u(2)^2)*(1+0.25*u(3)) + 0.1*u(3)*u(1)^3];

% Step 4: Solve the ODE.
for i=1:n-1
    k1 = h*f(u(:,i));
    k2 = h*f(u(:,i) + k1/2);
    k3 = h*f(u(:,i) + k2/2);
    k4 = h*f(u(:,i) + k3);
    u(:,i+1) = u(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
end
end