
a = 0;
b=2;
alpha =1;
beta=1;
e = 10^(-6);

h =0.2;
N= (b-a)/h;
x = a:h:b;

%The non-linear function
f = @(x,y) y.^2/(x+2);

%Intial condition
y0 = ones(1,N);

for i=1:N
    if i == 1
        r_con(i) = feval(f,x(i),y0(i))-alpha;
    elseif i == N
        r_con(i) = feval(f,x(i),y0(i))-beta;
    else
         r_con(i) = feval(f,x(i),y0(i));
    end
end
%The R vector
r_vec = r_con';

%Setting up matrix
upperdiag  = ones(1,N-1);
subdiag  = ones(1,N-1);
maindiag = -2*ones(1,N);

%Matrix A
A = diag(upperdiag,-1)+diag(maindiag)+diag(subdiag,1);

for i = 1:N
    
    Y1 = inv(A)*r_vec;
   
    %fprintf('Y%d = %.10f\n', i, Y1)
     if abs(Y1 - y0) <= e
         break
     end
     disp(i)
     y0 = Y1;
end





% 
% g =  @(x) 2*x +6 ;
% x0 = 1;
% e =10^(-3);
% n = 10;
% 
% for i = 1:n
%     x1 = g(x0);
%     fprintf('x%d = %.10f\n', i, x1)
%     if abs(x1 - x0) < e
%         break
%     end
%     x0 = x1;
% end