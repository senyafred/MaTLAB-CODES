% Function applying the modified Euler method to a given system of
% equations   

function [t,Y]=Frederick_HW7_modified_Euler(fun,tspan,h,y_init)

Y(:,1)=y_init;
t=tspan(1):h:tspan(2);

for n=1:length(t)-1
    Ybar = Y(:,n) + h*feval(fun,t(n),Y(:,n));
    Y(:,n+1) = .5*(Y(:,n) + Ybar + h*feval(fun,t(n+1),Ybar));
end

