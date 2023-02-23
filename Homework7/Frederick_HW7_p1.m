%%Solving y"+xy-3y=3x   y(0) = 1 y(2) = 5 
%using the Shooting method

% Set the parameters of the problem:
clear all;

%%Boundary condition
y0 = 1;  beta=5;

x0=0;                          %Initial value of x
xf=2;                          % final Value of x
h=0.1;                         % step size 

x=[x0 xf];                     % this is the vector of x values

yint=[1 0];                    %Initial conditions for 1st Auxillary
yint2=[0 1];                   % Initial conditions for 2nd Auxillary


% Solving Auxilliary IVPs using the Modified Euler method:
[x1,y1]=Frederick_HW7_modified_Euler(@Auxi1,x,h,yint);

[x2,y2]=Frederick_HW7_modified_Euler(@Auxi2,x,h,yint2);

%Finding parameter
theta = (beta - y1(1,end))/y2(1,end);

%The BVP solution
y= y1 + theta.*y2;     

%Plotting the Solution 
plot(x1,y(1,:),'b','linewidth',1.2)
xlabel('x')
ylabel('y')

function dy = Auxi1(x, y)
%First Auxilliary Equation
  dy = [y(2); -x*y(2)+3*y(1)+3*x];
end

function dy = Auxi2(x, y)
%Second Auxilliary Equation
  dy = [y(2); -x*y(2)+3*y(1)];
end