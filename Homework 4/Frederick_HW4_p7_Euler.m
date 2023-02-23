%%Euler Method
function [x,y] = Frederick_HW4_p7_Euler(F,yint,xspan,h)

y(1)=yint;
x=xspan(1):h:xspan(end);

for i=1:length(x)-1
    y(i+1)=y(i)+h*feval(F,x(i),y(i));
end
end