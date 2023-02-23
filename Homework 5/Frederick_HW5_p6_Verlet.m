% % This code integrates on [x0,xf] a model equation 
%Verlet

function [x,y1,y2] = Frederick_HW5_p6_Verlet(f1,f2,y1_int,y2_int,xspan,h)
%%Note y2 =v
y1(1)=y1_int;
y2(1)=y2_int;
x=xspan(1):h:xspan(end);

for i=1:length(x)-1
    y1(i+1) = y1(i) + h*y2(i) +(h^2/2)*f2(x(i),y1(i),y2(i));
    
    y2(i+1)= y2(i)+(h/2)*(f2(x(i),y1(i))+f2(x(i+1),y1(i+1),y2(i)));
     
end
end







