%%The Implicit Method

function [x,y1,y2] = Frederick_Test1_p4_Implicit(yint1,yint2,xspan,h)
%%Note y2 is the same as v

y1(1)=yint1;
y2(1)=yint2;
x=xspan(1):h:xspan(end);

for i=1:length(x)-1
    y1(i+1)=(y1(i)+h*y2(i))*(1/(1+h^2));   %%The method
    
    y2(i+1)=(y2(i)-h*y1(i))*(1/(1+h^2));  
end
end