%%Solving x^3y"'+xy'-y= -3 + Inx   y(1) = 1 y'(2) = 1/2  y"(2)=1/4 
%using the Shooting method

clear all;

%%Boundary condition
y0 = 1;  beta =1/2; gamma =1/4;

x0=1;                      %Initial value of x
xf=2;                      % final Value of x
h=0.02;                    % step size 

x=[x0 xf];                 % this is the vector of x values


yint=[1 0 0];              %Initial conditions for 1st Auxillary
yint2=[0 1 0];             % Initial conditions for 2nd Auxillary
yint3=[0 0 1];             % Initial conditions for 3rd Auxillary

% Solving Auxilliary IVPs using the Modified Euler method:
[x1,U]=Frederick_HW7_modified_Euler(@Auxi1,x,h,yint);

[x2,V]=Frederick_HW7_modified_Euler(@Auxi2,x,h,yint2);

[x3,W]=Frederick_HW7_modified_Euler(@Auxi3,x,h,yint3);

%Linear system
A = [V(2,end) W(2,end); V(3,end) W(3,end)];
B = [beta - U(2,end); gamma - U(3,end)];

%Solving the linear system
para = A\B;

%Findin the parameter
y= U + para(1).*V + para(2).*W;

%Plotting the Solution
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

function dy = Auxi3(x, y)
%Third Auxilliary Equation
  dy = [y(2); y(3); (1/(x^3))*(-x*y(2)+y(1))];
end