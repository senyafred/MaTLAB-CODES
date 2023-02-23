%% PROBLEM 1

%Coefficients of y
c1 = -1/10;                                    
c2 = (1/10*cos(sqrt(10)) - 1/10)/(sin(sqrt(10)));

%Coefficients of z
c11 = -1/9.9;                                    
c22 = (1/9.9*cos(sqrt(9.9)) - 1/9.9)/(sin(sqrt(9.9)));

%% The Analytical solutions for y and z
 y_analytical = @(x) c1*cos(sqrt(10)*x) + c2*sin(sqrt(10)*x) + 1/10;
 z_analytical = @(x) c11*cos(sqrt(9.9)*x) + c22*sin(sqrt(9.9)*x) + 1/9.9;

%Plot of the solution
subplot( 1, 2, 1);   fplot(y_analytical,[0 1])
xlabel('x')
ylabel('y')
title('Analytical Solution for y')  
subplot( 1, 2, 2 );  fplot(z_analytical,[0 1])
xlabel('x')
ylabel('z')
title('Analytical Solution for z') 



