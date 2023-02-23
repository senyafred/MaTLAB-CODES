
%% Problem 5 q3

%Plot of r versus (beta_h)

clear all;

h = 0.05;
rCoeff = 2;
r = 1/(rCoeff*h);
theta = 0.5;

beta_h = linspace(0,pi,100);

rho_exact1 = exp(-r*beta_h.^2);

rho1 = Rho(beta_h,r,theta);  %Function 1

rCoeff = 4;
r = 1/(rCoeff*h);
theta = 0.5;

rho_exact2 = exp(-r*beta_h.^2);

rho2 = Rho(beta_h,r,theta);  %Function 2

rCoeff = 2;
r = 1/(rCoeff*h);
theta = 1;

rho3 = Rho(beta_h,r,theta);  %Function 3

rCoeff = 4;
r = 1/(rCoeff*h);
theta = 1;

rho4 = Rho(beta_h,r,theta);  %Function 4

hold on
% plot(beta_h,rho1,beta_h,rho_exact1, beta_h,rho_exact2,beta_h,rho2,beta_h,rho3,beta_h,rho4); 
plot(beta_h,rho1,'LineWidth',1);
plot(beta_h,rho_exact1,"--",'LineWidth',1);
plot(beta_h,rho_exact2,"-.",'LineWidth',1);
plot(beta_h,rho2,"-o",'LineWidth',1);
plot(beta_h,rho3,":",'LineWidth',1);
plot(beta_h,rho4,".",'LineWidth',1);
xlabel("beta h")
ylabel("rho")
legend("1/2h, theta = 0.5","Exact1","Exact2", "1/4h, theta = 0.5","1/2h, theta = 1","1/4h, theta = 1")
hold off

function rho = Rho(beta_h,r,theta)

for i = 1:length(beta_h)
    
    z = 4*r*sin(beta_h(i)/2)^2;
    
    rho(i)  = (1 - (1-theta)*z)/(1 + theta*z);
end
end






