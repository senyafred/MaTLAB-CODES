%% MATLAB CODE FOR FTCS METHOD
% PROJECT WORK

clc;

xmin = 0;        %Minimum value of x
xmax = 1;        %Maximum value of x
Tmax= 0.5;       %Maximum value of t

kappa = 0.009;   %stepsize in t
h = 0.01;        % Stepsize in x

a=1;             % Wave speed

%Space and Time interval
x = xmin:h:xmax;
t = 0:kappa:Tmax;

M = length(x);
N = length(t);

%Courant-Friedrichs-Lewy
c = a*(h/kappa);   

% Gaussian Inital condition
u0 = exp(-200*(x-0.25).^2);
v0 = exp(-200*(x-0.25).^2);
u=u0;
v=v0;
unp1 = u0;
unp2 = v0;

for i =1:N-1
    
    for j =2:M-1
        
         unp1(j) = u(j) - c/2*(u(j+1) - u(j-1));
         unp2(j) = v(j) - c/2*(v(j+1) - v(j-1));
         
    end
    
    %update u
    u=unp1;
    v=unp2;
    
    %The Exact solution
    exact = exp(-200*(x - 0.25 - a*t(i)).^2);
    
    %Plot Solution
    plot(x,exact,'r-');        hold on;
    plot(x,u,'bo-',x,v,'k--');  hold off;
    axis([xmin xmax -0.5 1.5])
    xlabel('x','fontsize',14)
    ylabel('x','fontsize',14)
    title(sprintf('time = %1.3f',t(i)), 'fontsize',14)
    legend("Exact","Numerical Method")
    shg    %Show current figure
    
    % pause   

end

















% clc;
% 
% xmin = 0;        %Minimum value of x
% xmax = 1;        %Maximum value of x
% Tmax= 0.5;       %Maximum value of t
% 
% kappa = 0.009;   %stepsize in t
% h = 0.01;        % Stepsize in x
% 
% a=1;             % Wave speed
% 
% %Space and Time interval
% x = xmin:h:xmax;
% t = 0:kappa:Tmax;
% 
% %Courant-Friedrichs-Lewy
% c = a*(h/kappa);   
% 
% % Gaussian Inital condition
% u0 = exp(-200*(x-0.25).^2);
% u=u0;
% unp1 = u0;
% 
% 
% for i =1:length(t)-1
%     
%     for j =2:length(x)-1
%          unp1(j) = u(j) - c/2*u(j+1) + c/2*u(j-1);
%     end
%     
%     %update u
% 
%     u=unp1;
%     
%     %Calculate the exact solution
%     exact = exp(-200*(x - 0.25 - a*t(i)).^2);
%     
%     %Plot Solution
%     plot(x,exact,'r-');
%     hold on;
%     plot(x,u,'bo-','markerfacecolor','b');
%     hold off;
%     axis([xmin xmax -0.5 1.5])
%     xlabel('x','fontsize',14)
%     ylabel('x','fontsize',14)
%     title(sprintf('time = %1.3f',t(i)), 'fontsize',14)
%     legend("Exact","Numerical Method")
%     shg    %Show current figure
%     %pause(dt)
%     % pause   
% 
% end




% clc;
% 
% xmin =0;      %Minimum value of x
% xmax = 1;     %Maximum value of x
% N = 100;      %Number of nodes -1
% dt = 0.009;   %timestep
% t=0;          %time
% tmax= 0.5;    %Maximum value of t
% v = 1;        %Velocity
% 
% %discritise the domain
% dx = (xmax - xmin)/N;
% x = xmin-dx:dx:xmax+dx;
% 
% % Set initial condition
% uO = exp(-200*(x-0.25).^2);
% u=uO;
% unp1 = uO;
% 
% nsteps = tmax/dt;
% 
% for n =1: nsteps
%     
%     u(1) = u(3);
%     u(N+3) = u(N+1);
%     
%     for i =2:N+2
%         unp1(i) = u(i) -v*dt/dt*(u(i) -u(i-1));
%     end
%     
%     %update t and u
%     t=t+dt;
%     u=unp1;
%     
%     %Calculate the exact solution
%     exact = exp(-200*(x - 0.25 - v*t).^2);
%     
%     %Plot Solution
%     plot(x,exact,'r-');
%     hold on;
%     plot(x,u,'bo-','markerfacecolor','b');
%     hold off;
%     axis([xmin xmax -0.5 1.5])
%     xlabel('x','fontsize',14)
%     ylabel('x','fontsize',14)
%     title(sprintf('time = %1.3f',t), 'fontsize',14)
%     legend("Exact","Numerical Method")
%     shg
%     %pause(dt)
%      pause   
% 
% end
