% % This code integrates on [x0,xf] a model equation 
%Modified Euler

function [x,y1,y2] = Frederick_HW5_p6_modified_Euler(f1,f2,y1_int,y2_int,xspan,h)
%%Note y2=v

y1(1)=y1_int;
y2(1)=y2_int;
x=xspan(1):h:xspan(end);

halfh=h/2;

for i=1:length(x)-1
    z1 = y1(i) + h*f1(x(i),y1(i),y2(i));
    z2 = y2(i) + h*f2(x(i),y1(i),y2(i));
    y1(i+1)= y1(i)+halfh*(f1(x(i),y1(i),y2(i))+f1(x(i+1),z1,z2));
    y2(i+1)= y2(i)+halfh*(f2(x(i),y1(i),y2(i))+f2(x(i+1),z1,z2));  
end
end







