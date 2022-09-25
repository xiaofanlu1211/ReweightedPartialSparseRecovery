function [E,X,x,K] = WPSTRidge(A,D,b,lambda,gamma)

% ============================================================
% Goal: to solve the linear system Ax=b using the proposed
%       algorithm when partial support set D of x known a 
%       priori, where A' is computed using least-squares.
%
% Inputs:
%   A = matrix of size m*n
%   D = expect dense part
%   b = vector of size m*1
%   lambda = thresholding parameter
%   gamma = ridge parameter
%
% Outputs:
%   E = the value of objective function F(x)
%   X = the solution at each iteration
%   x = the approximate solution of linear system Ax=b
%   K = the number of inner iteration
% Authors: Xiaofan Lu, Huimei Ma, and Linan Zhang
% Date: Thursday, Sept 1, 2022
% Reference: Linan Zhang and Hayden Schaeffer,
% On the Convergence of the SINDy Algorithm
% https://epubs.siam.org/doi/abs/10.1137/18M1189828
% ============================================================

% parameters
[m,n] = size(A); err = 1; k = 2;     
% Maximum number of iterations 
MaxIt = n;                                                
% Expect sparse part
F = 1:n;                                      
S0 = setdiff(F,D);
% Auxiliary matrix
I = eye(m);
if isempty(D) == err
    P = I;
else
    P = I - A(:,D)*((A(:,D)'*A(:,D))\(A(:,D)'));  
end

% track the large indices and the solution in each step
LargeIndex = zeros(length(S0),MaxIt);        
SupportSet = zeros(length(S0),MaxIt);        
X = zeros(n,MaxIt);                           

% the initial step
x = zeros(n,1);
x(S0) = ((P*A(:,S0))'*(P*A(:,S0))+gamma*eye(length(S0)))\((P*A(:,S0))'*(P*b));                   
x(D) = A(:,D)\(b - A(:,S0)*x(S0));
X(:,1) = x;                             
% S = { j in S0: |x_j| >= lambda }
S = find( abs(x)>=lambda);
S = intersect(S,S0);
% sort S such that |x_S(j+1)| >= |x_S(j)|
[~,ind] = sort(abs(x(S)),'descend');
SupportSet(1:length(ind),1) = S;

% the iterative step
while k<=MaxIt && err
    
    % xnew = argmin ||PAy-Pb|| subject to supp(y) in S
    y = ((P*A(:,S))'*(P*A(:,S))+gamma*eye(length(S)))\((P*A(:,S))'*(P*b));                
    xnew = zeros(n,1);       
    xnew(S) = y;            
    
    % S = { j in S0 : |x_j| >= lambda }
    S = find( abs(xnew)>=lambda );
    S = intersect(S,S0);
    % sort S such that |x_S(j+1)| >= |x_S(j)|
    [~,ind] = sort(abs(xnew(S)),'descend');
    LargeIndex(1:length(ind),k) = S(ind);
    SupportSet(1:length(ind),k) = S;
    
    % check whether S^{k} = S^{k-1}
    err = ~isequal(SupportSet(:,k-1), SupportSet(:,k));
    
    % update
    x(S0) = xnew(S0);                     
    x(D) = A(:,D)\(b - A(:,S0)*x(S0));
    X(:,k) = x;
    k = k+1;
    
end
 
% inner iteration number
X(:,k:end) = [];
k = size(X,2);
K = k-1;
% find position of nonzeros of x
NX = X;
NX(NX~=0) = 1;

% norm of PA
normPA = norm(P*A(:,S0),2);

% the value of objective function F_wp
E = zeros(k,1);
for ii=1:k
    E(ii) = norm(P*A(:,S0)*X(S0,ii)-P*b,'fro')^2 / (normPA^2)+ gamma*norm(X(S0,ii),2)^2 + dot(lambda(S0).^2,NX(S0,ii));
end
end
