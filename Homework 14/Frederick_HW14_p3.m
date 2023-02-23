%%  Problem 3
clear;

Lmax = 5;     %Maximum length
Tmax = 3;   %Maximum time
h = 0.2;
kappa = h/4;  %Step size in t
gamma = 0.125;

%Space and time intervals
t = 0:kappa:Tmax;
x = -5:h:Lmax;

%x1 = -5+h:h:Lmax-h;

M = length(x);
N = length(t);

%Boundary condition
g_o = 0;
g_1 = 2;

%Initial condition
u0 = (4/pi)*atan(exp(7*x));

r = (kappa)/h^2;

%Starting method 12.12
for k=2:M-1
uhalf(k-1) = r/2*u0(k-1) + (1-2*r/2)*u0(k) + r/2*u0(k+1);
end
%  figure(3)
%   plot(x(2:M-1),uhalf,'b -.', 'LineWidth',1.7) 
%   pause
 
maindiag = (1 + (kappa*gamma)/h^2)*ones(1,M-2); 

% % %Vector
uold = u0(2:end-1);
%  
for i= 2:N
  %% Semi Implicit
    
    upperdiagonal =((-kappa*gamma)/(2*h^2) - kappa*(uhalf - 1)/(4*h)); 
    upperdiagonal = upperdiagonal(1:end-1);
    
    subdiagonal =((-kappa*gamma)/(2*h^2) + kappa*(uhalf - 1)/(4*h));
    subdiagonal = subdiagonal(2:end);
  
%     %Matrix B diagonals
     d = diag(-upperdiagonal,-1)+diag(2-maindiag)+diag(-subdiagonal,1);
     
     %vector b with boundary condtions
     b_vector = [2*((kappa*gamma)/(2*h^2) - kappa*(uhalf(1) - 1)/(4*h))*g_o,...
                 zeros(1,M-4),...
                 2*((kappa*gamma)/(2*h^2) + kappa*(uhalf(end) - 1)/(4*h))*g_1];
    
     %Final B matrix
      B = d*uold' + b_vector;
    
%     %Solving using Thomas
      unew = thomas(subdiagonal,maindiag,upperdiagonal,B);
       
      uold= unew;
     
%       plot(x(2:end-1), unew,'b -.', 'LineWidth',1.7) 
%       pause
     

end
% %Final solution
 u_final = [g_o,unew,g_1];

% %The exact solution
 for i= 1:M
    yexact(i) = 1 + tanh(x(i)/(2*gamma));     
 end
% 
% %Plot solution
figure(1)
 plot(x, u_final,'b -.', 'LineWidth',1.7) 
 xlabel('x')
 ylabel('y')
% 
 figure(2)
 plot(x,yexact,'b -.', 'LineWidth',1.7)
 xlabel('x')
 ylabel('y')
