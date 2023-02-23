%%  Problem 4

clear;

Lmax = 1;     %Maximum length
Tmax = 0.5;   %Maximum time

%Boundary condition
g_o = 0;
g_1 = 0;

h = 0.1;
kappa = input('enter the value for parameter kappa:  ');  %Step size in t

%Space and time intervals
t = 0:kappa:Tmax;
x = 0:h:Lmax;

M = length(x);
N = length(t);

%Initial condition
u0 = sin(pi*x);

r = kappa/h^2;

%Vectors
subdiagonal = (-r/2)*ones(1,M-3);
upperdiagonal = (-r/2)*ones(1,M-3);
maindiag = 1-(r/2*(-2))*ones(1,M-2);

for i= 1:N-1
    
   uold = u0;
   
    for j = 2:M-1
        
     d(j-1)= (r/2)*uold(j+1) + (1 + (r/2)*(-2))*uold(j) + (r/2)*uold(j-1);

    end
    
      unew = thomas(subdiagonal,maindiag,upperdiagonal,d);
   
      unew = [g_o, unew ,g_1];
      
        u0 = unew;
end

%The exact solution
for i= 1:N
    for j = 1:M
    yexact(i,j) = sin(pi*x(j))*exp(-(pi^2)*t(i));   
    end   
end

hold on;
plot(x, unew,'b -.', 'LineWidth',1.7) 
xlabel('x')
ylabel('y')
legend(" k = 0.1","k = 0.05","k = 0.025","k = 0.01")
hold off;

% hold on;
% plot(x,yexact(end,:)- unew,'b -.', 'LineWidth',1.7)
% xlabel('x')
% ylabel('y')
% legend("error @ k = 0.1","error @ k = 0.05"," error @ k = 0.025","error @ k = 0.01","Exact")
% hold off;
