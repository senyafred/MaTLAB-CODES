%%Solving y" = y^2/(2+x) y(0) = 1 y(2) = 1 with Tol = 10^-3
%using the Shooting method

clear all;
% Set the parameters of the problem:

x0=0;                          %Initial value of x
xf=2;                          % final Value of x
h=0.01;                        % step size 

x=[x0 xf];                     % this is the vector of x values
y0=1;                          %initial condition 
        
theta1= -15;                   %initial guess for theta (Deep Solution)
theta2= -10;   

%theta1= 0;                   %initial guess for theta (Shallow Solution)
%theta2= -5;  

iteration = 1;

while(abs(theta2-theta1) > 10^(-3))               %Setting Tolerance
    
    %Solving Auxilliary IVPs using the Modified Euler method:
    [x1,yt1]=Frederick_HW7_modified_Euler(@Auxi,x,h,[y0 theta1]);      
    [x2,yt2]=Frederick_HW7_modified_Euler(@Auxi,x,h,[y0 theta2]);  
    
    y1 = yt1(1,end) -1;             %| F(theta_k) |
    y2 = yt2(1,end) -1;
    
    %Secant Method
    m = theta2 - ((theta2-theta1)*(y2)/((y2)-(y1)));  
    theta1 = theta2;
    theta2 = m;
    
    iteration = iteration + 1;
    
end

% Show number of iteration
disp(iteration);

[xt,yt1]=Frederick_HW7_modified_Euler(@Auxi,x,h,[y0 m]);

%Plotting the solution
figure
plot(xt,yt1(1,:),'b','linewidth',1.2); 
xlabel('x')
ylabel('y')
title('Deep solution when theta1 =-10, theta2 = -15')

function dy=Auxi(x,y) 
%BVP function
    dy = [y(2);(1/(2+x))*y(1)^2];
end

