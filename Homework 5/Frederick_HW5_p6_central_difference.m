% % This code integrates on [x0,xf] a model equation 
%Central Difference


function [x,y,v] = Frederick_HW5_p6_central_difference(f,y1_int,y2_int,xspan,h)

%% Note y2 = v

y(1)=y1_int;
v(1)=y2_int;
secondcondi = y2_int;

x=xspan;
y(2)=y(1) + h*secondcondi + (h^2/2)*f(x,y(1));

for i=2:length(x)-1
   
    y(i+1)= 2*y(i) - y(i-1) + h^2*feval(f,x(i),y(i));   %The method
     
    v(i) = (y(i+1) - y(i-1))/(2*h);
    
    if i == length(x)-1
       v(i+1) = (y(i+1)-y(i))/h + h/2*(-y(i+1));
    end
   
end
end








