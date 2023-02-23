%%The Implicit Method
function [x,y] = Frederick_HW4_p7_Implicit(yint,xspan,h)
k0 = -20;       % Constant from the equation
y(1)=yint;
x=xspan(1):h:xspan(end);

for i=1:length(x)-1
    y(i+1)=(1/(1-k0*h))*y(i);    
end
end
