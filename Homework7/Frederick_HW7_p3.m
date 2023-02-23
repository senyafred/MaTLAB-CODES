%%Solving x^3y"'+xy'-y= -3 + Inx   y(1) = 1 y'(2) = 1/2  y"(2)=1/4 
%using the Shooting method

clear all;
%%Boundary condition
y0 = 1;  beta =1/2; gamma =1/4;

x0=2;                            %Initial value of x
xf=1;                            % final Value of x
h=-0.02;                         % step size 

x=[x0 xf];                       % this is the vector of x values

yint=[0 beta gamma];             %Initial conditions for 1st Auxillary
yint2=[1 0 0];                   %Initial conditions for 2nd Auxillary
        
% Solving Auxilliary IVPs using the Modified Euler method:
[x1,U]=Frederick_HW7_modified_Euler(@Auxi1,x,h,yint);

[x2,V]=Frederick_HW7_modified_Euler(@Auxi2,x,h,yint2);

%Find parameter
theta = (y0 - U(1,end))/V(1,end);

%The BVP solution
y= U + theta.*V;  

%Plot the solution
plot(x1,y(1,:),'r', x1,y(2,:),'k--',x1,y(3,:),'bo','linewidth',1.1)
xlabel('x')
ylabel('y')
legend(' Y1', 'Y2','Y3');

function dy = Auxi1(x, y)
%First Auxilliary Equation
  dy = [y(2); y(3); (1/(x^3))*(-x*y(2)+y(1)+(-3+log(x)))];
end

function dy = Auxi2(x, y)
%Second Auxilliary Equation
  dy = [y(2); y(3); (1/(x^3))*(-x*y(2)+y(1))];
end
