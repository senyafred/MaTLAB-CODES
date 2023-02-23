% Usage: [y t] = abm4(f,a,b,ya,n) or y = abm4(f,a,b,ya,n)
% Adams-Bashforth-Moulton 4-th order predictor-corrector method for initial value problems
% It uses 
% Adams-Bashforth 4-step method as a precdictor,
% Adams-Moulton 3-step method as a corrector, and
% Runge-Kutta method of order 4 as a starter
%
% Input:
% f - Matlab inline function f(t,y)
% a,b - interval
% ya - initial condition
% n - number of subintervals (panels)
%
% Output:
% y - computed solution
% t - time steps
%
% Examples:
% [y t]=abm4(@myfunc,0,1,1,10);         % here 'myfunc' is a user-defined function in M-file
 %y=abm4(inline('sin(y*t)','t','y'),0,1,1,10);
 
 f = @(y) sin(y);
 [y,t]=abm4(f,0,1,1,10); 


function [y, t] = abm4(f,a,b,y0,n)
h = (b - a) / n;

y(1) = y0;
t(1) = a;

m = min(2,n);

for i = 1 : m % start-up phase, using Runge-Kutta of order 4
    s(i) = f(y(i));
    y(i+1) = y(i) + h*(s(i));
    
end

for i = m + 1 : n % P-C method
    s(i) = f(y(i));
    y(i+1) = y(i) + h/2*(3 * s(i) - s(i-1)); %  Predictor

    y(i+1) = y(i) + h/2*(s(i) + f(y(i+1))); %% Corrector
end
end