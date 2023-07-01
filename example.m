clear all
% Construct matrix A, vector x and b
m = 10;
n = 5;
eta = randn(m,1)*1;                % random noise~N(0,1)
A = randi([-10 10],m,n);           % random matrix
% A = normc(A);
x = ones(n,1);
x(1:2) = 10;                       % true x
b = A*x + eta;                     
epsilon = 0.01;                    % parameters
gamma = 0;

%% Example
% clc

% A = [-6 -7 1 -8 0;
    9 -4 1 5 6;
    0 -5 -9 -3 8;
    0 -4 -3 -5 4;
    1 9 -3 -5 0;
    8 -2 7 7 -2;
    0 -2 -4 1 0;
    -9 -1 -3 -6 1;
    6 -4 6 -5 7;
    0 8 2 8 0];
% b = [-135.99, 62.58, -54.87, -44.93, 93.72, 71.93, -22.70, -105.99, 28.51, 88.82]';

Ds = {1, 1:2, 1:3, 1:4};          
for i=1:length(Ds)
    D = Ds{i};                      
    lambda0 = [8*ones(length(D),1);0.8*ones(n-length(D),1)];
    [E,X,x,K] = WPSTRidge(A,D,b,lambda0,gamma);
    disp(x)
    [x,lambda,Xout,Xin,kout,kin,F] = RWPSTRidge(A,D,b,lambda0,epsilon,gamma);
    disp(x)
end

%%
% save('example1.mat','A','x','b','eta');
