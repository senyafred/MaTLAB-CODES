clear all;

%% Problem 4
%Solving Tridiagonal Matrix Ay = r

M = 1000;          %Matrix dimension
%M = 5000;

%Setting up matrix A
upperdiag  = -1*ones(1,M-1);
subdiag  = -1*ones(1,M-1);
maindiag = 2*ones(1,M);

%Matrix A (regular matrix)
A = diag(upperdiag,-1)+diag(maindiag)+diag(subdiag,1);

%Vector r
r = (-ones(1,M)).^(0:M-1)';

%Using Thomas algorithm to solve the system
for k = 1 : 200
    tic
    Y_thomas = thomas(subdiag ,A, upperdiag ,r);
    tt(k) = toc;
end
t_thomas = sum(tt)    %Computational time


%% Using Matlab Solver for the system
%Case I
for j = 1:200
    tic
    Y_full = A \ r;  %Solving the system
    t(j) = toc;
end
tfull = sum(t)    %Computational time


%Case II
%Matrix A (Sparse matrix)
a  = -1*ones(M,1);
c  = -1*ones(M,1);
b = 2*ones(M,1);

%Vector r
r1 = (-ones(1,M)).^(0:M-1);

sparseA = spdiags(a,-1,M,M)+spdiags(b,0,M,M)+spdiags(c,1,M,M);

for i=1:200
    tic
    Y_sparse = sparseA/r1;      %Solving the system
    T(i)= toc; 
end
tsparse = sum(T)    %Computational time




