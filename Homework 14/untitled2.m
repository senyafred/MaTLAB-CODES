%This is code for pr1 hw 14

clear all

%parameteres 

h = 0.02;
k = h/2;
x = 0 :h: 1; 
M = length(x);
t = 0 :k: 2;
N = length(t); 
p = ones(N);
r = k/h^2;

%Boundary value conditions
g= 0*ones(1,N);
q(1:N)= pi;
p = ones(1,N);

%at t=0

U(1:M)= sin(pi*x(1:M));

a = (-r/2)*ones(1,M-2);
c(1)=-r;
c(2:M-2) = (-r/2);
b(1)= (1 + r*(1-h*p(2)));
b(2:M-1) = (1+r);
 
for n= 2:N
    
   Uold(1:M) = U(1:M);

    R(1) = (r)*Uold(2) + (1 - r*(1 - h*p(1)))*Uold(1) - r*h*(q(1) + q(2));
    
    for m = 2:M-1
        R(m)= (r/2)*Uold(m+1) + (1 - r)*Uold(m) + (r/2)*Uold(m-1);
    end
    R(M-1) = (r/2)*Uold(M-2) + (1 - r)*Uold(M-1) + (g(1) + g(2))*r/2;
   
     U=thomas(a,b,c,R);
     U=[U,g(n)];

end


figure(1)
plot(x,U)
title('plot of the error when t=2')
xlabel('x')
ylabel('U')
