%%Solving y"+xy-3y=3x   y(0) = 1 y(2) = 5 
%using the Shooting method

% Set the parameters of the problem:
clear all;

%%Boundary condition
alpha = -15/19;  beta=0;

x0=2;                          %Initial value of x
xf=3;                          % final Value of x
h=0.1;                         % step size 

x=[x0 xf];                     % this is the vector of x values

yint=[0 -15/19];                    %Initial conditions for 1st Auxillary
yint2=[1 0];                   % Initial conditions for 2nd Auxillary


% Solving Auxilliary IVPs using the Modified Euler method:
[x1,y1]=Frederick_HW7_modified_Euler(@Auxi1,x,h,yint);

[x2,y2]=Frederick_HW7_modified_Euler(@Auxi2,x,h,yint2);
 
 %Finding parameter
 theta = (beta - y1(1,end))/y2(1,end);
 
 %The BVP solution
 y= y1 + theta*y2;     
 
 
 % The exact solution
yexact = @(x) (114*x-19*x^2 - 5*x^3 - 36)/(38*x);
for k = 1:length(x1)
    yexact1(k) = yexact(x1(k));
end
 
 %Plotting the Solution 
 plot(x1,y(1,:),'k',x1,yexact1,'r','linewidth',1.2)
% xlabel('x')
% ylabel('y')

function dy = Auxi1(x, u)
%First Auxilliary Equation
  dy = [u(2); (2/x^2)*u(1)+((x-6)/x^2)];
end

function dy = Auxi2(x, u)
%Second Auxilliary Equation
  dy = [u(2); (2/x^2)*u(1)];
end





