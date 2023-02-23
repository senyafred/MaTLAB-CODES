% Homework 15 - Problem 9:

clear all;

K = 0.002;
h = 0.1;
tmax=0.5;
ymax = 2.5;

r = K/h^2;

x=0:h:1;
y=0:h:ymax;
t=0:K:tmax;

M=length(x);   
L=length(y);
    
N=length(t);   
    
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
        
    %Boundary Condition    
         U(L,:)= x(1:M)*exp(-(2*pi/ymax)^2*t(n));
         U(1,:)= x(1:M)*exp(-(2*pi/ymax)^2*t(n));
         U(:,M)= cos(2*pi*y(1:L)./ymax)*exp(-(2*pi/ymax)^2*t(n));
         U(:,1)= sin(2*pi*y(1:L)./ymax)*exp(-(2*pi/ymax)^2*t(n));  
         
         
         G0 = 1/2*(Uold(2:L-1,1) + U(2:L-1,1)) + r/4*((Uold(3:L,1) - 2*Uold(2:L-1,1))...
             + Uold(1:L-2,1)) - (U(3:L,1) - 2*U(2:L-1,1) + U(1:L-2,1));
                 
         G1 = 1/2*(Uold(2:L-1,M) + U(2:L-1,M)) + r/4*((Uold(3:L,M) - 2*Uold(2:L-1,M))...
             + Uold(1:L-2,M)) - (U(3:L,M) - 2*U(2:L-1,M) + U(1:L-2,M));
         
         Ustar(2:L-1,1) = G0;
         Ustar(2:L-1,M) = G1;
        
         
     for ell = 2: L-1       
             B = (Uold(ell,2:M-1) + (r/2)*(Uold(ell+1,2:M-1) -2*Uold(ell,2:M-1) + Uold(ell-1,2:M-1))...
                 + r/2*[G0(ell-1),zeros(ell-1),G1(ell-1)])';

     end
              
    % Tridiagonal matrix in x-direction
     subdiagonal = (-r/2)*ones(1,M-3);
     upperdiagonal = (-r/2)*ones(1,M-3);
     maindiag = (1+r)*ones(1,M-2);
%     
%     % Thomas algorithm to solve intermetiate solution 
%           
      u_star = thomas(subdiagonal,maindiag,upperdiagonal,B);
%   
     for m = 2: M-1       
             C = (u_star(2:L-1,m) + (r/2)*(u_star(2:L-1,m+1) -2*u_star(2:L-1,m) + u_star(2:L-1,m-1))...
                 + r/2*[G0(ell-1),zeros(ell-1),G1(ell-1)])';
     end

    % Tridiagonal matrix in y-direction
    subdiagonal = (-r/2)*ones(1,M-3);
    upperdiagonal = (-r/2)*ones(1,M-3);
    maindiag = (1+r)*ones(1,M-2);


    unew = thomas(subdiagonal,maindiag,upperdiagonal,C);
    
    U = Unew;
%      
%         
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




