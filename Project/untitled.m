clc;
clear;

xmin =0;      %Minimum value of x
xmax = 1;     %Maximum value of x
Tmax = 0.05;
dx = 0.01;
dt = 0.009;


x = xmin:dx:xmax;
t = 0:dt:Tmax;

M =length(x);
N =length(t);

a=1;    %Velocity

c = a*(dx/dt);   %Courant

% Set initial condition
uO = exp(-200*(x-0.25).^2);
u=uO;
unp1 = uO;


for n =1:length(t)-1
    
    
    for i =2:length(M)
        
        %unp1(i) = u(i) -v*dt/dt*(u(i) -u(i-1));
        unp1(i) = u(i) - c/2*u(i+1) + c/2*u(i-1);
    end
    
%     update t and u
%     t=t+dt;
      u=unp1;
%     
%     Calculate the exact solution
       exact(n,i) = exp(-200*(x(i) - 0.25 - a*t(n)).^2);
%     
%     Plot Solution
%     plot(x,exact,'r-');
%     hold on;
%     plot(x,u,'bo-','markerfacecolor','b');
%     hold off;
%     axis([xmin xmax -0.5 1.5])
%     xlabel('x','fontsize',14)
%     ylabel('x','fontsize',14)
%     title(sprintf('time = %1.3f',t), 'fontsize',14)
%     legend("Exact","Numerical Method")
%     shg    %Show current figure
%     pause(dt)
%     pause   

end
