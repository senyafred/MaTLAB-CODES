
clear all;

a = 0;
b=2;
alpha =1;       %BC at a
beta=1;         %BC at b
e = 10^(-6);    %Tolerance
c=0.75;

h = 0.1;       %Step size
N = (b-a)/h;   %Number of iteration
x = a:h:b;     % X interval

%Setting up matrix A
upperdiag  = ones(1,N-2);
subdiag  = ones(1,N-2);
maindiag = -2*ones(1,N-1);

%Matrix A
A1 = diag(upperdiag,-1)+diag(maindiag)+diag(subdiag,1);
I = c*eye(N-1)*h^2;
A = A1 - I;

%The non-linear function
f = @(x,y) y.^2/(x+2);

%Intial condition
y0 = ones(N-1,1);

for k = 1:N-1
    
    for j = 1:N-1
        if j ==1
            r(j,1)= h^2*f(x(j),y0(j))-alpha;
        elseif j == N-1
            r(j,1)= h^2*f(x(j),y0(j))-beta;
        else
            r(j,1)= h^2*f(x(j),y0(j));
        end           
    end
     
     y11 = (A - c*I*h^2)\(r - c*h^2*y0);
      
      if norm(y11 - y0) <= e
          break
      end
      
     y0 = y11;
     disp(k)
end

y = [alpha;y11;beta];

%Plot solution
plot(x,y)














% 
% clear all;
% 
% a = 0;
% b=2;
% alpha =1;       %BC at a
% beta=1;         %BC at b
% e = 10^(-6);    %Tolerance
% 
% h = 0.1;       %Step size
% N = (b-a)/h;   %Number of iteration
% x = a:h:b;     % X interval
% 
% %Setting up matrix A
% upperdiag  = ones(1,N-2);
% subdiag  = ones(1,N-2);
% maindiag = -2*ones(1,N-1);
% 
% %Matrix A
% A = diag(upperdiag,-1)+diag(maindiag)+diag(subdiag,1);
% 
% %The non-linear function
% f = @(x,y) y.^2/(x+2);
% 
% %Intial condition
% y0 = [alpha;ones(N-1,1);beta];
% 
% for k = 1:N-1
%     
%     for j = 1:N-1
%         if j ==1
%             r(j,1)= h^2*f(x(j),y0(j))-alpha;
%         elseif j == N-1
%             r(j,1)= h^2*f(x(j),y0(j))-beta;
%         else
%             r(j,1)= h^2*f(x(j),y0(j));
%         end           
%     end
%     
%       y11 = A\r;
%       
%       y1 = [alpha;y11;beta];
%       
%       if norm(y1 - y0) < e
%          break
%       end
%       
%      y0 = y1;
% end
% %y1 = [alpha;y11;beta];
% 
% %Plot solution
% plot(x,y0)






