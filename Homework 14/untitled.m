%%  Problem 2
clear;

Lmax = 5;     %Maximum length
Tmax = 3;   %Maximum time
h = 0.02;
kappa = h/4;  %Step size in t
gamma = 0.125;

%Space and time intervals
t = 0:kappa:Tmax;
x = -5:h:Lmax;

x1 = -5+h:h:Lmax-h;

M = length(x1);
N = length(t);

% x_halfmin = (x(1:M-2) + x(2:M-1))/2;
% x_halfmax = (x(2:M-1) + x(3:M))/2;

%Boundary condition
g_o = 0;
g_1 = 2;

%Initial condition
u0 = (4/pi)*atan(exp(7*x1));

r = (kappa)/h^2;

% u0 = u0(2:end-1);

%using Simple Euler 12.12
for m=2:M-1
u_n_onehalf(m-1) = r/2*u0(m-1) + (1-2*r/2)*u0(m) + r/2*u0(m+1);
end
 
% %Vectors
  upperdiagonal =((-kappa*gamma)/(2*h^2) - kappa*(u_n_onehalf - 1)/4*h);  
  maindiag(2:M-1) = (1 + (kappa*gamma)/h^2); 
  subdiagonal =((-kappa*gamma)/(2*h^2) + kappa*(u_n_onehalf - 1)/4*h);
  
uold = u0(2:end-1);
 
for i= 2:N
  %% Crank Nicholsen
    
    %Matrix A (regular matrix)
    d = diag(-upperdiagonal,-1)+diag(2-maindiag)+diag(-subdiagonal,1);

    B = d*uold';
    

%      unew = thomas(subdiagonal,maindiag,upperdiagonal,B);
% %       
%       u0 = [g_o,unew,g_1];'

end