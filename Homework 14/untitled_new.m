%This is the code for problem 3 HW 14

clear all

%parameteres 

h = 0.1;   % have to decide on that using hw 13             
k = h/4;    %for now I just use the same values as in pr 2
x= -5 :h:5;
M = length(x);
t = 0 :k: 3;
N = length(t); 
r = k/h^2;
gamma = 0.125;

%Boundary value conditions
g1= 0*ones(1,N);
g2= 2*ones(1,N);

% u at t=0
U(1:M-2) = (4/pi)*atan(exp(7*x(2:M-1)));
U= [g1(1),U,g2(1)];

%calculate u at level 0+1/2 in time to start method 14.50
%using Simple Euler 12.12

for m=2:M-1
u_n_onehalf(m) = r/2* U(m-1) + (1-2*r/2)*U(m) + r/2*U(m+1);
end
      

for n = 2:N
    uold(1:M-2)= U(2:M-1);

    %superdiags and the main diag for matrices B and A respectively
       c_aux(1:M-2) = (k)*( gamma/(2*h^2) + (u_n_onehalf(1:M-2) -1)/(4*h));
       a_aux(1:M-2) = (k)*( gamma/(2*h^2) + (u_n_onehalf(1:M-2) -1)/(4*h)); 
       c(1:M-3)= -c_aux(1:M-3); 
       a(1:M-3)= -a_aux(2:M-2);
       b_aux= 1-(k)*(gamma/(2*h^2))*ones(1,M-2);
       b = 1 + (k)*(gamma/(2*h^2))*ones(1,M-2);
 
       B = spdiags(a_aux',-1,M-2,M-2) + spdiags(b_aux',0,M-2,M-2) + spdiags(c_aux',1,M-2,M-2);
       beta = [0*ones(1,M-3),(k)*( gamma/(2*h^2) + (g2(n) -1)/(4*h))*g2(n)];
       R= B*uold' + beta;
       U= thomas(a,b,c,R);
       u_n_onehalf = (3/2)*U - uold/2;
       U= [g1(n),U,g2(n)];    
       
%        figure(777)
%        plot (x,U)
%        pause

end



u_exact(1:M) = 1 + tanh(x(1:M)./2*gamma);
ue(1:M)=sech(gamma*x(1:M)).* tanh(gamma*x(1:M));

figure(1)
plot (x,U)
title('plot of the numerical sol')

figure (2)
plot (x, u_exact-U)
title('error ')

figure (3)
plot (x, ue)


