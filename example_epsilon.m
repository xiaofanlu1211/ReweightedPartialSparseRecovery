% =========================================================================
% Goal: Compare the influence of epsilon on RWPSTRidge Algorithm
% 
% Authors: Xiaofan Lu, Huimei Ma, and Linan Zhang
% Date: Sept, 2022
% =========================================================================
% Construct matrix A and vector b

m = 256; n = 100;                   % size of matrix A
D = 1:5;                            % dense part/prior known support set of x
eta = 0.01;                         % noise level
lambda0 = 0.02*ones(n,1);           % initial threshold parameter
gamma = 0;                          % ridge parameter
P = zeros(20,5);                    
x0 = zeros(n,1);
for k = 1:20                        % sparsity of x
    S = 1:k;                        % position of nonzeros of x
    x0(S) = 1;                      % construct ''true'' vector x
    % the iteration step
    for i = 1:100                   
        A = randn(m,n);             % matrix A with i.i.d. Gaussian entries
        b = A*x0 + eta*randn(m,1);  % construct b with noise                         
        for q = -3:1                
            epsilon = 10^q;         % 5 different epsilon values
            % find x
            x = RWPSTRidge(A,D,b,lambda0,epsilon,gamma);
            % check whether supp(x) = supp(x0)\cup D
            suppx = find(x~=0);
            err = isequal(suppx,union(S,D)');
            % Record the number of perfect recoveries.
            if err == 1 
                P(k,q+4) = P(k,q+4)+1; 
            end
        end
    end
end
% Calculate recovery probability
P = P/100; 