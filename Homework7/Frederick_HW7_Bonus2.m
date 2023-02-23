%%Solving y" = 30^2(y -1 + 2x) y(0) = 1 y(1.62) = -2.24 
%using Multiple Shooting methods

clear all;
%%Boundary condition
y0 = 1;  b =1.62;  beta=-2.24;

x0=0;                          %Initial value of x
xf=b/2;                        % final Value of x
h=0.01;                        % step size 

x=[x0 xf];                   % this is the vector of x values

yint=[1 0];                  %Initial conditions for 1st Auxillary
yint2=[0 1];           

% Solving Auxilliary IVPs in the First Interval
% using the Modified Euler method:
[x1,U1]=Frederick_HW7_modified_Euler(@Auxi1,x,h,yint);
[x2,V1]=Frederick_HW7_modified_Euler(@Auxi2,x,h,yint2);

%%  Solving Auxilliary IVPs in the Second Interval
xx=[b/2 b];   %this is the vector of x values

[xx1,U2] =Frederick_HW7_modified_Euler(@Auxi1,xx,h,[0 0]);
[xx1,V21]=Frederick_HW7_modified_Euler(@Auxi2,xx,h,[1 0]);
[xx1,V22]=Frederick_HW7_modified_Euler(@Auxi2,xx,h,[0 1]);
%%
%Linear system
A = [V1(1,end) -1 0; V1(2,end) -0  -1; 0 V21(1,end) V22(1,end)];
B = [-U1(1,end); -U1(2,end);beta - U2(1,end)];

%Solving linear system
para = A\B;

%Two solutions to the BVP
w1= U1 + para(1).*V1; 
w2= U2 + para(2).*V21 + para(3).*V22;

%Creating a vector for the two solutions
W_vec = [w1 w2];
x_vec = [x1 xx1];

% The exact solution
yexact = @(x)(1-2*x);

%Plotting Exact and Numerical Solution
plot(x_vec,W_vec(1,:),'b',x_vec,yexact(x_vec),'r--','linewidth',1.3); 
xlabel('x')
ylabel('y')
legend(' Numerical', 'Exact');

function dy = Auxi1(x, y)
%First Auxilliary Equation
  dy = [y(2); 30^2*(y(1) - 1 + 2*x)];
end

function dy = Auxi2(x, y)
%Second Auxilliary Equation
  dy = [y(2); 30^2*y(1)];
end