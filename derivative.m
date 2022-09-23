function dXdt = derivative(Xtilde,h)
% ============================================================
% Description:
%   The function outputs numerical derivative
%
% Inputs:
%   Xtilde = given data x
%   h = time step
%
% Outputs:
%   dXdt = 1-th time derivative of Xtilde
% Authors: Xiaofan Lu, Huimei Ma, and Linan Zhang
% ============================================================

% Calculate numerical derivative
[m,n] = size(Xtilde);                                        %calculate Xtilde size          
dXdt = zeros(m,n);                  
dXdt(1,:) = (Xtilde(2,:) - Xtilde(1,:))./h;            
dXdt(2:m-1,:) = (Xtilde(3:m,:) - Xtilde(1:m-2,:))./(2*h);    %central difference method
dXdt(m,:) = (Xtilde(m,:) - Xtilde(m-1,:))./h;
end