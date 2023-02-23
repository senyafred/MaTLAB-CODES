%% MATLAB CODE FOR LAX METHOD
% PROJECT WORK

clc;

xmin = -1;        %Minimum value of x
xmax = 2;        %Maximum value of x
Tmax= 1;       %Maximum value of t

% kappa = 0.01;   %stepsize in t
% h = 0.01;        % Stepsize in x
h = 0.05;
kappa = 0.9*h;   %stepsize in t
        % Stepsize in x

a=1;             % Wave speed

%Space and Time interval
x = xmin:h:xmax;
t = 0:kappa:Tmax;

M = length(x);
N = length(t);

%Courant-Friedrichs-Lewy
c = a*(kappa/h); 
%c = 1; 

% Gaussian Inital condition
u0 = exp(-200*(x-0.25).^2);
v0 = exp(-200*(x-0.25).^2);
u=u0;
v=v0;
unp1 = u0;
unp2 = v0;

for i =1:N-1
    
    for j =2:M-1
        
         unp1(j) = 1/2*(u(j+1) + u(j-1)) - c/2*(u(j+1) - u(j-1));
         unp2(j) = 1/2*(v(j+1) + v(j-1)) + c/2*(v(j+1) - v(j-1));
         
    end
    
    %update u
    u=unp1;
    v=unp2;
    
    %The Exact solution
    exact1 = exp(-200*(x - 0.25 - a*t(i)).^2);
    exact2 = exp(-200*(x - 0.25 + a*t(i)).^2);
    
    %Plot Solution
     plot(x,exact1,'r-',x,exact2,'r-');    hold on;
    plot(x,u,'b-',x,v,'k-');  hold off;
    axis([xmin xmax -0.5 1.5])
    xlabel('x','fontsize',17)
    ylabel('x','fontsize',17)
    title(sprintf('time = %1.3f',t(i)), 'fontsize',14)
    legend("Exact1","Exact2","Numerical Method","Numerical Method2")
    shg    %Show current figure
    
    % pause   

end
%Error
figure(5000)
plot(x,exact1-u,x,exact2-v,'r-');


