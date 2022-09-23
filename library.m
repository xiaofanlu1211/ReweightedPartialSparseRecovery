function [D,A,L] = library(u,p,s)

% =========================================================================
% Goal: construct dictionary matrix for find underlying equations
% 
% Inputs:
% u = data used to construct dictionary matrix A and D
% p = the maximum degree of polynomials in D
% s = order of known polynomial

% Outputs:
% A = a library of polynomials of degree at most p
% D = index of known polynomials of order s 
% L = cell that serves as the legend of A
% Authors: Xiaofan Lu, Huimei Ma, and Linan Zhang
% Date: June, 2022
% =========================================================================
[m,n]=size(u);         % returns the dimension of u
x1=[ones(m,1) u];      % Add a column vector of all 1s to u;(size:m*n+1)
l1 = cell(n+1,1);      % array with (n+1)*1
if n==2                
    l1{2} = 'u';       % let the second line of l1 be u
    l1{3} = 'v';       % let the third line of l1 be v
else
    for ii=1:n
        l1{ii+1} = strcat('u',num2str(ii));% let the ii+1 line be uii
    end
end
C = 1:(n+1);
for ii=1:p-1
    C=combvec(C,1:(n+1));     % Create all combinations of 1:n+1;(size:p*(n+1)^p)
end
C=sort(C,1);                  %arrange the column elements in C in ascending order
C=unique(C','rows');          %delete duplicate columns in C and return as rows;(size:c*p)
nA=size(C,1);                 %returns the length of the first(row) dimension of C;(size:c*1)
A=ones(m,nA);                 %Create a 1 matrix;(size:m*c)
L=cell(nA,1);                 %create a array;(size:c*1)
for ii=1:nA                   
    for jj=1:p                
        A(:,ii)=A(:,ii).*x1(:,C(ii,jj));      %Multiply the corresponding columns.
        L{ii}=strcat(L{ii},l1{C(ii,jj)});     %paste the elements in l1 into L.
    end
end
D = sum(C~=1,2);              %find the position of polynomials of order s
D = find(D==s);         
L{1} = num2str(1);              %make the first line be 1