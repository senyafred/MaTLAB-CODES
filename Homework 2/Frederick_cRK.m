%%Classical Runge-Kutta Method

function [x,Y] = Frederick_cRK(F,yint,xspan,h)

Y(1)=yint;
x=xspan(1):h:xspan(2);

for i=1:(length(x)-1)                              % calculation loop   
    k_1 = feval(F,x(i),Y(i));
    k_2 = feval(F,x(i)+(h/2),Y(i)+(h/2)*k_1);
    k_3 = feval(F,x(i)+(h/2),Y(i)+(h/2)*k_2);      %Classical Runge-Kutta
    k_4 = feval(F,x(i)+h,Y(i)+k_3*h);

    Y(i+1) = Y(i) + (h/6)*(k_1+2*k_2+2*k_3+k_4);  % main equation
end

  