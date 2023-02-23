% setup and  parameters
Nx = 50;  % number of mesh points
dx = 10/Nx; % spactial mesh 
Tend =1;% solve from t=0..tmax

c = 0.2*ones(Nx-1,1); %  vector  of ones , size Nx-1
cfl=0.95;
dt = cfl*dx/max(c);% time step
R = round(Tend/dt); %  number of timesteps
x = linspace(dx,10-dx,Nx-1);% create spatial coordinates vector for plots

% set up differentiation matrices 
e = ones(Nx-1,1); % vector of ones of same size x
Dx  =(spdiags([-e e],[-1 1],Nx-1,Nx-1)); % 1st order matrix
Dxx = (spdiags([e -2*e e], [-1 0 1], Nx-1,Nx-1)); % 2nd order matrix

% initialization and initial conditions,   u = zero at boundaries 
u = exp(-100 * (x-5).^2)';

for n = 1:R % iteratre, step in time
   u(:,n+1) = u(:,n) - 0.5*dt/dx * c.*(Dx*u(:,n)) + 0.5*(dt/dx)^2 * c.^2.*(Dxx*u(:,n)) ...
      + 0.125*(dt/dx)^2 * c.*(Dx*c).*(Dx*u(:,n));
end

(plot(x,u(:,R+1),'o - r'));xlabel('x'); ylabel('u'); title('Solution @ t=1')