%%Solving y" = -y y(0)=0  y'(0) = 1 using different methods

% Set the parameters of the problem:
clear all;

f = @(x,y) -y;                  % function


f1 = @(x,y1,y2) y2;              % System of ODE
f2 = @(x,y1,y2) -y1; 

x0=0;                          %Initial value of x
xf=1000;                       % final Value of x
h=0.2;                         % step size for second order
h1=0.5;                        % Step size for fourth order

x=[x0:h:xf];                   % this is the vector of x values


y1(1)=0;                      %Initial conditions
y2(1)=1;                       %% Where y2 is V

% Solve the equation using the Verlet method:
[y_Vt,y_V1,v_V1]=Frederick_HW5_p6_Verlet(f1,f2,y1(1),y2(1),x,h);
% Hamiltonian Error Verlet method:
H_verlet = 0.5 -1/2*(y_V1.^2 + v_V1.^2);
% Solve the equation using the Modified Euler method:
[mt,Y1_mE,Y2_mE]=Frederick_HW5_p6_modified_Euler(f1,f2,y1(1),y2(1),x,h);
% Hamiltonian Error Modified Euler
H_modified = 0.5 - 1/2*(Y1_mE.^2 + Y2_mE.^2);
% Solve the equation using the cRK method:
[rkt,Y1_cRK,Y2_cRK]=Frederick_HW5_p6_cRK(f1,f2,y1(1),y2(1),x,h1);
%Hamiltonian Error cRK
H_cRK = 0.5 - 1/2*(Y1_cRK.^2 + Y2_cRK.^2);
% Solve the equation using the Central Difference:
[ct,y_CD,v_CD]=Frederick_HW5_p6_central_difference(f,y1(1),y2(1),x,h);
%Hamiltonian Error Central Difference
H_central = 0.5 - 1/2.*(y_CD.^2 + v_CD.^2);
% Solve the equation using ODE45:
options=odeset('AbsTol',0.002);
[t1_ode45,y1_ode45]=ode45(@Frederick_HW5_p6_vectorF,[x0,xf],[0;1],options);
%Hamiltonian Error ODE45
H_ode1 = 0.5 - 1/2*(y1_ode45(:,1).^2 + y1_ode45(:,2).^2);
% Solve the equation using ODE45:
options=odeset('AbsTol',0.003);
[t2_ode45,y2_ode45]=ode45(@Frederick_HW5_p6_vectorF,[x0,xf],[0;1],options);
%Hamiltonian Error ODE45
H_ode2 = 0.5 -1/2*(y2_ode45(:,1).^2 + y2_ode45(:,2).^2);


%%%%%%%%%%%%%%
%Phase Plane
%%%%%%%%%%%%%%

figure(50601);
% Plot the Verlet solution:  ---------->
subplot(3,2,1);  plot(y_V1, v_V1);
axis([ -1.1*max(abs(y_V1)) 1.1*max(abs(y_V1)) ...
       -1.1*max(abs(v_V1)) 1.1*max(abs(v_V1)) ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('Verlet-1','fontsize',14)
% Plot the modified Euler solution:  ---------->
subplot(3,2,2);  plot(Y1_mE, Y2_mE);
axis([ -1.1*max(abs(Y1_mE)) 1.1*max(abs(Y1_mE)) ...
       -1.1*max(abs(Y2_mE)) 1.1*max(abs(Y2_mE)) ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('modified Euler','fontsize',14)
% Plot the cRK solution:  ---------->
subplot(3,2,3);  plot(Y1_cRK, Y2_cRK);
axis([ -1.1*max(abs(Y1_cRK)) 1.1*max(abs(Y1_cRK)) ...
       -1.1*max(abs(Y2_cRK)) 1.1*max(abs(Y2_cRK)) ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('classical Runge--Kutta','fontsize',14)
% Plot the central-difference solution:  ---------->
subplot(3,2,4);  plot(y_CD, v_CD);
axis([ -1.1*max(abs(y_CD)) 1.1*max(abs(y_CD)) ...
       -1.1*max(abs(v_CD)) 1.1*max(abs(v_CD)) ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('central-difference','fontsize',14)
% Plot the ode45 solution with AbsTol = 0.002:  ---------->
subplot(3,2,5);  plot(y1_ode45(:,1), y1_ode45(:,2));
axis([ -1.01 1.01 -1.01 1.01 ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('ode45 with Tol=0.002','fontsize',14)
% Plot the ode45 solution with AbsTol = 0.003:  ---------->
subplot(3,2,6);  plot(y2_ode45(:,1), y2_ode45(:,2));
axis([ -1.01 1.01 -1.01 1.01 ])
axis('equal')
xlabel('y', 'fontsize', 12);   ylabel('v', 'fontsize', 12);
title('ode45 with Tol=0.003','fontsize',14)

%%%%%%%%%%%%%%
%Hamiltonians
%%%%%%%%%%%%%%

figure(50602);
% Plot the Verlet solution:  ---------->
subplot(3,2,1);  plot(y_Vt, H_verlet);
% axis([ -1.1*max(abs(y_V1)) 1.1*max(abs(y_V1)) ...
%        -1.1*max(abs(v_V1)) 1.1*max(abs(v_V1)) ])
% axis('equal')
xlabel('time', 'fontsize', 12);   ylabel('Error (H)', 'fontsize', 12);
title('Hamiltonian Verlet-1','fontsize',14)
% Plot the modified Euler solution:  ---------->
subplot(3,2,2);  plot(x, H_modified);
% axis([ -1.1*max(abs(Y1_mE)) 1.1*max(abs(Y1_mE)) ...
%        -1.1*max(abs(Y2_mE)) 1.1*max(abs(Y2_mE)) ])
% axis('equal')
xlabel('time', 'fontsize', 12);   ylabel('Error (H)', 'fontsize', 12);
title('Hamiltonian Modified','fontsize',14)
% Plot the cRK solution:  ---------->
subplot(3,2,3);  plot(rkt, H_cRK);
% axis([ -1.1*max(abs(Y1_cRK)) 1.1*max(abs(Y1_cRK)) ...
%        -1.1*max(abs(Y2_cRK)) 1.1*max(abs(Y2_cRK)) ])
% axis('equal')
xlabel('time', 'fontsize', 12);   ylabel('Error (H)', 'fontsize', 12);
title('Hamiltonian cRK','fontsize',14)
% Plot the central-difference solution:  ---------->
subplot(3,2,4);  plot(ct, H_central);
% axis([ -1.1*max(abs(y_CD)) 1.1*max(abs(y_CD)) ...
%        -1.1*max(abs(v_CD)) 1.1*max(abs(v_CD)) ])
% axis('equal')
xlabel('time', 'fontsize', 12);   ylabel('Error (H)', 'fontsize', 12);
title('Hamiltonian Central Difference','fontsize',14)
% Plot the ode45 solution with AbsTol = 0.002:  ---------->
%subplot(3,2,5);  plot(Y_ode45_A(:,1), Y_ode45_A(:,2));
subplot(3,2,5);  plot(t1_ode45, H_ode1);
% axis([ -1.01 1.01 -1.01 1.01 ])
% axis('equal')
xlabel('time', 'fontsize', 12);   ylabel('Error (H)', 'fontsize', 12);
title('Hamiltonian Ode1','fontsize',14)
% Plot the ode45 solution with AbsTol = 0.003:  ---------->
%subplot(3,2,6);  plot(Y_ode45_B(:,1), Y_ode45_B(:,2));
subplot(3,2,6);  plot(t2_ode45, H_ode2);
% axis([ -1.01 1.01 -1.01 1.01 ])
% axis('equal')
xlabel('time', 'fontsize', 12);   ylabel('Error (H)', 'fontsize', 12);
title('Hamiltonian Ode2','fontsize',14)



