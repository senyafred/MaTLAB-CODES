% % This code integrates on [x0,xf] a model equation 
%cRK

function [x,y1,y2] = Frederick_HW5_p6_cRK(f1,f2,y1_int,y2_int,xspan,h)
%%Note y2=v

y1(1)=y1_int;
y2(1)=y2_int;
x=xspan(1):h:xspan(end);


for i=1:(length(x)-1)                              % calculation loop   
    k_11 = feval(f1,x(i),y1(i),y2(i));
    k_12 = feval(f2,x(i),y1(i),y2(i));
    k_21 = feval(f1,x(i)+(h/2),y1(i)+(h/2)*k_11,y2(i)+(h/2)*k_12);
    k_22 = feval(f2,x(i)+(h/2),y1(i)+(h/2)*k_11,y2(i)+(h/2)*k_12);
    k_31 = feval(f1,x(i)+(h/2),y1(i)+(h/2)*k_21,y2(i)+(h/2)*k_22); 
    k_32 = feval(f2,x(i)+(h/2),y1(i)+(h/2)*k_21,y2(i)+(h/2)*k_22);
    k_41 = feval(f1,x(i)+h,y1(i)+k_31*h,y2(i)+k_32*h);
    k_42 = feval(f2,x(i)+h,y1(i)+k_31*h,y2(i)+k_32*h);

    y1(i+1) = y1(i) + (h/6)*(k_11+2*k_21+2*k_31+k_41);  % main equation
    y2(i+1) = y2(i) + (h/6)*(k_12+2*k_22+2*k_32+k_42);
end
end








