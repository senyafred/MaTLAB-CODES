% This code integrates on [x0,xf] a model equation 
%   y' = sqrt(y)
% with initial condition y(x0)=y0.
% The integration is done by the Modified Euler method.


% Set the parameters of the problem:
clear all;

f1 = @(x,y1,y2) y2;              % function
f2 = @(x,y1,y2) -y1; 


x0=0;                          %Initial value of x
xf=1000;                          % beginning and ending values of x
h=0.2;                         % step size in x
                                 % initial value of 
w=1;                                
x=[x0:h:xf];                   % this is the vector of x values


%y1(1)=0;
%y2(1)=1;

% for i=1:length(x)
%     y1(i)= exp(((h*w^2)/2)*x(i))*abs(sin(x(i)));
%     y2(i)= exp(((h*w^2)/2)*x(i))*abs(cos(x(i)));
%     
% end


y1= exp(((h*w^2)/2)*xf)*abs(sin(xf));
y2= exp(((h*w^2)/2)*xf)*abs(cos(xf));
  
%
H = 1/2*(y2^2 + y1^2)
%plot(y2,y1)

Diff =  H - 0.5