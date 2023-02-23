% %%  Problem 2
clear;

Lmax = 1;     %Maximum length
Tmax = 0.5;   %Maximum time
h = 0.02;
kappa = h/8;  %Step size in t

%Space and time intervals
t = 0:kappa:Tmax;
x = 0:h:Lmax;

M = length(x);
N = length(t);

x_halfmin = (x(1:M-2) + x(2:M-1))/2;
x_halfmax = (x(2:M-1) + x(3:M))/2;

%Boundary condition
g_o = 0;
g_1 = 0;

%Initial condition
u0 = sin(pi*x);

r = kappa/h^2;

%Vectors
subdiagonal = -r/2.*(1 + x_halfmin);
subdiagonal = subdiagonal(2:end);

upperdiagonal = -r/2.*(1 + x_halfmax);
upperdiagonal = upperdiagonal(1:end-1);

uold = u0(2:end-1);

for i= 2:N
  %% Crank Nicholsen
   
   for j = 2:M-1
       
       maindiag(j-1) = 1 + r/2*(1+x_halfmax(j-1)) + r/2*(1 + x_halfmin(j-1)) + kappa*x(j)^2/(2*(1 + t(i)^2));
            
   end
   
    %Matrix B
    d = diag(-upperdiagonal,-1)+diag(2-maindiag)+diag(-subdiagonal,1);
   
    B = d*uold';

    unew = thomas(subdiagonal,maindiag,upperdiagonal,B);
%     plot(x(2:end-1), unew,'b -.', 'LineWidth',1.7) 
%     pause
      
    uold= unew;
      
end
%Final Solution
u_final = [g_o,unew,g_1];

plot(x, u_final,'b -.', 'LineWidth',1.7) 
xlabel('t')
ylabel('U')
