%%Solving y" = 30^2(y -1 + 2x) y(0) = 1 y(1.62) = -2.24 
%using the Shooting methods

clear all;

%%Boundary condition
y0 = 1;  beta=-2.24;

x0=0;                 %Initial value of x
xf=1.62;              % final Value of x
h=0.01;               % step size 

x=[x0 xf];            % this is the vector of x values
          
yint=[1 0];           %Initial conditions for Auxillary IVP
yint2=[0 1];

% Solving Auxilliary IVPs using the Modified Euler method:
[x1,y1]=Frederick_HW7_modified_Euler(@Auxi1,x,h,yint);

[x2,y2]=Frederick_HW7_modified_Euler(@Auxi2,x,h,yint2);

%Finding parameter
theta = (beta - y1(1,end))/y2(1,end);

%BVP solution
y= y1 + theta.*y2;    

% The exact solution
yexact = @(x)(1-2*x);

%Plotting Exact and Numerical Solution
plot(x1,yexact(x2),'r--',x1,y(1,:),'b','linewidth',1.3);   
xlabel('x')
ylabel('y')
legend(' Exact', 'Numerical');


function dy = Auxi1(x, y)
%First Auxilliary Equation
  dy = [y(2); 30^2*(y(1) - 1 + 2*x)];
end

function dy = Auxi2(x, y)
%Second Auxilliary Equation
  dy = [y(2); 30^2*y(1)];
end