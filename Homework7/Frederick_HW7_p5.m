%%Solving y" + (2sech^2x -lambda^2)y 

x0=-10;                          %Initial value of x
xf=10;                           % final Value of x
h=0.01;                          % step size for second order
                                
R=10; 
c=1;
x=[x0 xf];

dlamda = 0.01;                  %Step in lamda

%Intial conditon
yint = exp(-c*R);

for lamda = 0.1:dlamda:10
     yint_prime = lamda*yint;       %initial condtion
     
     %Using ODE45
    options=odeset('AbsTol',0.001);
    [t,G_lamda] = ode45(@(x,y) vectorF(x,y,lamda),x,[yint yint_prime],options);
    
    if (G_lamda(end,1)) > 0
        break
    end
    
    lamda_star = lamda;
end

%Final shooting
yint_prime = lamda_star*yint;
[x1,G_lamda] = ode45(@(x,y) vectorF(x,y,lamda_star),x,[yint yint_prime]);

%Plot the graph
 plot(x1,G_lamda(:,1))
    
function z = vectorF(x, y,lamda)
%BVP Function
    z(1) = y(2); 
    z(2) = -(2*(sech(x)).^2 - lamda^2)*y(1);
    
    z = transpose(z);
end