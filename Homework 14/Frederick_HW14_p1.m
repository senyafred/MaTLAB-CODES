%%  Problem 1

clear;

Lmax = 1;     %Maximum length
Tmax = 2;   %Maximum time
h = 0.02;
kappa = h/2;  %Step size in t

%Space and time intervals
t = 0:kappa:Tmax;
x = 0:h:Lmax;

M = length(x);
N = length(t);

%Boundary
p= ones(1,N);
q =pi*ones(1,N);
g1 =zeros(1,N);


%Initial condition
u0 = sin(pi*x);

r = kappa/h^2;

%Vectors
subdiagonal = (-r/2)*ones(1,M-2);

upperdiagonal = [-r, (-r/2)*ones(1,M-3)];

maindiag = [(1 +r*(1-h*p(2))), (1+r)*ones(1,M-2)];

for i= 1:N-1
  %% Crank Nicholsen
   uold = u0;
   
   for j = 2:M-1
        
      d(1)=  (1-r*(1-h*p(1)))*uold(1) + (r)*uold(2) +(-r*h*(q(1)+q(2))); 
                
        if j == M-1
            d(j)= (r/2)*uold(j-1) + (1-r)*uold(j) + ((r/2)*(g1(1)+g1(2)));
            
        else
             d(j)= (r/2)*uold(j+1) + (1-r)*uold(j) + (r/2)*uold(j-1);
        end
    end
    
      unew = thomas(subdiagonal,maindiag,upperdiagonal,d);
      
      u0 = unew;
end

%Final solution
unew = [unew ,g1(N)];

plot(x, unew,'b -.', 'LineWidth',1.7) 
xlabel('t')
ylabel('U')
