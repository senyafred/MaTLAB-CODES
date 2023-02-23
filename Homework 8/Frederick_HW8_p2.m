%% Problem 2

%Solving a linear system AX = B

A = [ -2  2  0  0  0;
       1 -2  1  0  0;
       0  2 -2  0  0;
       0  0  3 -2 -1
       0  0  0  4 -2];
   
% B matrix
b = [2;0;-2;-4;4];

%Solving the system
x = A \ b;

