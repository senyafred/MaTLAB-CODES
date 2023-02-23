% Homework 15 - Problem 2:

clear all;

K = 0.002;
h = 0.1;
tmax=0.5;
ymax = 2.5;

r = K/h^2;

x=0:h:1;
y=0:h:ymax;
t=0:K:tmax;

M=length(x);   % THIS M = [(M+1) OF THE NOTES], i.e. M_notes = thisM - 1.
L=length(y);

    
N=length(t);   % THIS N = [(N+1) OF THE NOTES], i.e. N_notes = thisN - 1.
    
% % Initial condition:
for m = 1:M
    for ell = 1:L
        U0(ell,m) = 10*sin(pi*x(m)).*sin(pi*y(ell)/ymax) + sin(2*pi*y(ell)/ymax).*(1-x(m))...
            + cos(2*pi*y(ell)/ymax).*x(m); 
    end
end

U = U0;
    % Advance the solution from the n-th to the (n+1)-st level
  for n = 1 : N-1
      
        Uold = U;
        
        for ell = 2: L-1
            for m = 2 : M-1                
                Unew(ell,m) = r*Uold(ell,m+1) + (1-4*r)*Uold(ell,m) + r*Uold(ell,m-1) ...
                    + r*Uold(ell+1,m) + r*Uold(ell-1,m);
            end
        end
        
         U(2:L-1,2:M-1) = Unew(2:L-1,2:M-1);
        
         U(L,:)= x(1:M)*exp(-(2*pi/ymax)^2*t(n));
         U(1,:)= x(1:M)*exp(-(2*pi/ymax)^2*t(n));
         U(:,M)= cos(2*pi*y(1:L)./ymax)*exp(-(2*pi/ymax)^2*t(n));
         U(:,1)= sin(2*pi*y(1:L)./ymax)*exp(-(2*pi/ymax)^2*t(n));
        
%         [Xi,Yi] = meshgrid(x,y);
%         surf(Xi,Yi,U)
%         pause;
        
 end
    
subplot( 1, 2, 1);
[Xi,Yi] = meshgrid(x,y);
surf(Xi,Yi,U)
title("Solution to U_t = U_xx + U_yy")

%The exact solution
yexact = @(t,x,y) 10*sin(pi*x).*sin((pi*y)/ymax).*exp(-(1+(1/ymax^2))*pi^2*t) + (sin((2*pi*y)/ymax).*(1-x)...
            + cos((2*pi*y)/ymax).*x).*exp(-((2*pi)/ymax)^2*t); 

subplot( 1, 2, 2);
[Xi,Yi] = meshgrid(x,y);
surf(Xi,Yi,yexact(t(end),Xi,Yi))
title("The Exact Solution")

figure(3)
surf(Xi,Yi,U - yexact(t(end),Xi,Yi))
title("Error")




