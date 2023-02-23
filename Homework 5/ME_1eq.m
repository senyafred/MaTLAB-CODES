% Function for Modified Euler Method for 1 equation:
% y'=fun(x,y1,y2)
%
% the user must input fun, initial value for y,
% a 2-dimensional vector for xspan, and the desired step size.


function [x,Y]=ME_1eq(fun,yint,xspan,h);

Y(1)=yint;
x=xspan(1):h:xspan(2);

for n=1:length(x)-1
    Ybar=Y(n)+h*feval(fun,x(n),Y(n));
    Y(n+1)=1/2*( Y(n)+h*feval(fun,x(n+1),Ybar) + Ybar);
end