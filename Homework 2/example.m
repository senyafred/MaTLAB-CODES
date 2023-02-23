%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sahebeh Dadboud: 
% Explicit Euler methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




y0 = 1;                  % Initial Condition
h = 0.2;                 % step size
t = 0:h:10;               % t goes from 0 to 10 
yexact = -t + 1/100;      % Exact solution 

ystar(1) = y0;           % Initial condition 
for i=1:(length(t)-1)
    f = (1-ystar(i)/100)*ystar(i);  % y'(t) = (1-y(t)/100)*y(t)
    ystar(i+1) = ystar(i) + f*h; 
    error_explic(i) = abs(ystar(i)-yexact(i));
end
plot(t,yexact,t,ystar,error_explic);
legend('Exact','Approximate','Error');
disp(error_explic);
