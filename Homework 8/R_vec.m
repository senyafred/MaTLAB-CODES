
a = 0;
b=2;
alpha =1;
beta=1;
e = 10^(-6);


h =0.1;
N= (b-a)/h;
x = a:h:b;

%The non-linear function
f = @(x,y) y.^2/(x+2);

%Intial condition
y0 =ones(1,N-1);

for i=1:N-1
    if i == 1
        r_con(i) = h^2*feval(f,x(i),y0(i))-alpha;
    elseif i == N-1
        r_con(i) = h^2*feval(f,x(i),y0(i))-beta;
    else
         r_con(i) = h^2*feval(f,x(i),y0(i));
    end
end
%The R vector
r_vec = r_con';

%Setting up matrix
upperdiag  = ones(1,N-2);
subdiag  = ones(1,N-2);
maindiag = -2*ones(1,N-1);

%Matrix A
A = diag(upperdiag,-1)+diag(maindiag)+diag(subdiag,1);

for i = 1:N
    
    Y1 = A\r_vec;
   
     if norm(Y1 - y0) <= e
         break
     end
     disp(i)
     y0 = Y1;
end
y = [alpha;y0;beta];

plot(x,y)
