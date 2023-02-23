
clear all;
clc;
close all;
%%%  DEFINING PARAMETERS GIVEN %%%
a = 0.65;       % alpha in ft^2/hr
w = 1;          % width of bar (feet)
dx = 0.05;      % step size in x-direction (feet)
nx = 1+w/dx;    % no of steps in x-direction
L = 1;          % length of bar (feet)
dy = 0.05;      % step size in y-direction (feet)
ny = 1+L/dy;    % no of steps in y-direction
tf = 0.1;       % Final time
dt = 0.002;     % Time Step
Nt = tf/dt+1;   % Total number of time steps
%%% DEFINING ADDITIONAL CONSTANTS TO SIMPLIFY EQUATION
Dx = a*dt/(dx^2); % Diffusion number(dx)
A = 2*(1+Dx);
Dy = a*dt/(dy^2); % Diffusion number (dy)
B = 2*(1-Dy);
C = 2*(1+Dy);
D = 2*(1-Dx);
% Define 'x'
for i=1:nx
    x(i)=(i-1)*dx;
end
% Define 'y'
for j=1:ny
    y(j)=(j-1)*dy;
end
% Define 't'
for k=1:Nt
    t(k)=(k-1)*dt;
end
%%% INITIAL CONDITIONS %%%
for k = 1:Nt
    for i = 2:nx
        for j = 2:ny
            T(i,j,k) = 0;
        end
    end
end
%%% BOUNDARY CONDITIONS %%%
for k=1:Nt
    for i = 1:nx
        for j = 1:ny
            if j==1
                T(i,j,k)    =  200;   % Lower wall
            elseif i==nx
                T(i,j,k)   =  0;    % Right Wall
            elseif i==1
                T(i,j,k)    =  200;   % Left Wall
            elseif j==ny
                T(i,j,k)   =  0;    % Top Wall
            end
        end
    end
end
%   DEFINING MATRIX ON RHS
U=zeros(nx-2,1);
Ty=zeros(ny-2,1);
for k=1:Nt-1
    for j=2:ny-1
        for i=2:nx-1
            if i==1
                U(i,j) = Dy*u0(i,j+1) + B*u0(i,j) + Dy*u0(i,j-1) + Dx*u0(i,j+1);
            elseif i==nx-1
                U(i,j) = Dy*u0(i,j+1) + B*u0(i,j) + Dy*u0(i,j-1) + Dx*u0(i,j+1);
            else
                U(i,j) = dy*u0(i,j+1) + u0(i,j) + dy*u0(i,j-1);
            end
        end
    end
    
    % Tridiagonal matrix in x-sweep
    % First we transform scalar constant B and Dx into column vectors
    BB_M = B*ones(nx-1,1);
    DDx_U = Dx*ones(nx-2,1);
    DDx_L = DDx_U;
    
    ma_x = diag(BB_M) + diag(-DDx_U,1) + diag(-DDx_L,-1);
    % Find the solution at the interior nodes 
    T_x = (inv(ma_x))*U;
    
    % Tridiagonal matrix in y-sweep
    CC_M = C*ones(ny-1,1);
    DDy_U = Dy*ones(ny-2,1);
    DDy_L = DDy_U;
    ma_y = diag(CC_M) + diag(-DDy_U,1) + diag(-DDy_L,-1);
    
    % to solve y sweep
    T_y = inv(ma_y)*T_x;
    
  
    % To start over
    T=T_y;
end
pcolor(T)


    %Vectors
%     subdiagonal = (-r/2)*ones(1,M-3);
%     upperdiagonal = (-r/2)*ones(1,M-3);
%     maindiag = (1+r)*ones(1,M-2);
