%%Solving y" = -y y(0)=0  y'(0) = 1 using different methods

% Set the parameters of the problem:
clear all;

%%Boundary condition
f1= @(x,y1,y2) y2; 
f2= @(x,y1,y2) -x*y2 + 3*y1 + 3*x;

f3=@(x,y1,y2) y2;
f4=@(x,y1,y2) -x*y2 + 3*y1;

y0 = 1;  beta=5;

x0=0;                          %Initial value of x
xf=2;                       % final Value of x
h=0.1;                         % step size for second order

%t=[x0:h:xf];                   % this is the vector of x values
x=[x0 xf];

% y1(1)=1; y2(1)=0;             %Initial conditions for 1st Auxillary
yint=[1 0];
% z1(1)=0; z2(1)=1;             %% Initial conditions for 2nd Auxillary
yint2=[0 1];

y1(1) =1;
y2(1) = 0;

% Solve the equation using the Modified Euler method:
[x1,y1,y2]=Frederick_HW5_p6_modified_Euler(f1,f2,1,0,x,h);
%u=z(:,1)

[x2,y3,y4]=Frederick_HW5_p6_modified_Euler(f3,f4,0,1,x,h);

 z = (beta-y1(end))/y3(end);
  

 %plot(x1,w1,'r')

 

