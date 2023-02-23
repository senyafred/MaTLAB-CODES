% homework 13 - Problem 5:
 clear all;
% analyzing the effect of smoothness of the initial data on
%  the order of the error in the simple explicit method for the Heat equation
%   u_t = u_xx  with the boundary conditions  u(0)=u(1)=0.

% Setting up the color table for the lines (needed for plotting below):
line_color(1,:)=[1 0 0];   % red
line_color(2,:)=[0 0 0];   % black
line_color(3,:)=[0 0 1];   % blue
line_style=char('-','--','-.');

theta=1;          % 0.5 => Crank-Nicolson method
                    % 1   => implicit Euler method
                    
p=1;                % parameter characterizing the order of the lowest derivative that
                    %  is discontinuous in the initial condition
                    
coeff_r=input(' r = 1/(coeff_r*h); thus, enter 2 or 4: coeff_r = ');   % kappa = r*h^2     

tmax=0.5;

disp(' ')
disp(['Value of p = ' int2str(p)])


for kk=1:3          % the loop changing the step size in  x
	h=0.05/2^(kk-1);
    r=1/(coeff_r*h);          %  K/h^2
	K=r*h^2;
	
	x=0:h:1;
	M=length(x);   % for historical reasons, "this M" = "M in the notes" + 1
	t=0:K:tmax;
	N=length(t);
    
    %Boundary conditions
    g_o = 0;
    g_1 = 0;
    
    % Initial condition:
	Mo4=round((M-1)/4);  % auxiliary number = M/4 (see below)
    if p == 0
	    U0 = [zeros(1,Mo4)  ones(1,M-2*Mo4)  zeros(1,Mo4)];
    else 
        U0 = [zeros(1,Mo4)  (cos(2*pi*(x(Mo4+1:M-Mo4)-1/2))).^p   zeros(1,Mo4)];
    end
    U = U0;
	
    % -----------------------------------------------------------------------------------------
 
    %Vectors
    subdiagonal = -(r*theta)*ones(1,M-3);
    upperdiagonal = -(r*theta)*ones(1,M-3);
    maindiag = (1+r*2*theta)*ones(1,M-2);

    for n= 1:N-1
         uold = U;
        for m = 2:M-1
            d(m-1)= (r-r*theta)*uold(m+1) + (1 + (r-r*theta)*-2)*uold(m) + (r-r*theta)*uold(m-1);
        end
    
         unew = thomas(subdiagonal,maindiag,upperdiagonal,d);
   
         Unew = [g_o, unew ,g_1];
         
         U = Unew;
    end

    
    %  ----------------------------------------------------------------------------------------
    
    % Save the coarsest x-grid:
    if kk==1
        xcoarse=x;
    end
    
    % Record the final solution for each  h:
    Uend(kk,:) = U(1:2^(kk-1):end);
    % It is not used for plotting, but for the calculation of gamma.
    
    % Plot the solution at t=tmax for each step size  h.
    %  First, prepare an auxiliary index  nr  for subplots.
    %  With this index, the code can be run for two different values of  r  (and for
    %  given  p and  theta), and the results can be plotted all in one figure.
    if coeff_r == 2
        nr=1;
    elseif coeff_r == 4
        nr=2;
    else
        nr = 4444;
        disp(' ')
        disp('You entered an unspecified value of  coeff_r !  <--------------')
        disp('Rerun the code with a specified value.    <--------------')
        disp(' ')
        break
    end
    if nr <= 2
        figure(130501); 
		subplot(2,1,nr);
    end        
	plot(x,U,'color',line_color(kk,:),'linestyle',line_style(kk,:),'linewidth',2)
    ylabel(['r=1/(' num2str(coeff_r) '*h), p=' num2str(p)], 'fontsize', 16)
    if nr == 2
        xlabel('X', 'fontsize', 12)
    elseif nr == 1
        title(['solution at T_{max}=' num2str(tmax) '; p=' num2str(p)], 'fontsize', 14)
    end
    
	hold on
    
    % This is the same plot as in figure 1 (except for the initial condition),
    %  but it is not a subplot, and it is not saved when we run different values of  r.
    
    figure(130502);
    
    plot(x,U,'color',line_color(kk,:),'linestyle',line_style(kk,:),'linewidth',2)
    xlabel('X', 'fontsize', 12);
    ylabel(['U(t=end)  for  r=1/(' num2str(coeff_r) ' h)'], 'fontsize', 14)
    title('Same as fig. 130501, but stand-alone','fontsize',14)
    hold on
    
    % Clear the arrays whose dimensions will be changing with  h
    %  except for the last values of h (they will be needed for plotting
    %  of the initial condition):
    if kk < 3-eps
        clear('x','t','U','Uu','R','a','b','c');
    end   
        
end
% Plot the initial condition scaled to the amplitude of the final solution:
% 
figure(130501); subplot(2,1,nr);
plot(x,U0*max(abs(U))/max(abs(U0)),'m:','linewidth',3)
              % One needs to scale the amplitude of the initial condition so that it 
              %  could be plotted along with the solution at t=tmax. The above scaling 
              %  assumes that max(U0)=1, as it should be.
legend('h', 'h/2', 'h/4', '@t=0, h/4, scaled')
hold off


% Compute the reference (i.e., very close to exact) solution 
% by the Explicit Euler with a very small Kappa:
%  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
KREF = h^2/10;  % use the latest saved h, which is the finest 
              % Automatically, the x-grid used here is the finest one.
rREF = KREF/h^2;
NREF = round(tmax/KREF)+1;
% Boundary confitions:
g0REF = zeros(NREF,1);
g1REF = zeros(NREF,1);
% Advance the solution from the n-th to the (n+1)-st level
%  by the Explicit Euler method.

% As the initial condition, use only the values of U inside the spatial boundaries:
UREF = U0;
%
% Compute the solution inside the boundaries at a new time level and
%  then augment the resulting matrix to include the solution at the boundaries:
for n = 1 : NREF-1
    UREFold = UREF;
    for m = 2 : M-1                
        UREFnew(m) = rREF*UREFold(m+1) + (1-2*rREF)*UREFold(m) + rREF*UREFold(m-1);
    end
    UREF = UREFnew;
    UREF(1) = g0REF(n+1);
    UREF(M) = g1REF(n+1);
end
%  ----------------------------------------------------------------------------------------
    
figure(130502); 
plot(x,UREF,'g--', 'linewidth',2)
legend('h', 'h/2', 'h/4', '"exact"')
hold off

% Compute gamma_h  (in the notations of the HW problem):
error_exponent=log((Uend(1,:)-Uend(2,:))./(Uend(2,:)-Uend(3,:)))/log(2);

% Plot the real and imaginary parts of gamma versus x:
figure(130503);
subplot(2,1,1), plot(xcoarse,real(error_exponent))
xlabel('X', 'fontsize', 14)
ylabel('Re (\gamma)', 'fontsize', 14)
subplot(2,1,2), plot(xcoarse,imag(error_exponent))
xlabel('X', 'fontsize', 14)
ylabel('Im (\gamma)', 'fontsize', 14)

maxRegamma____maxImgamma = ...
[max(real(error_exponent)), max(imag(error_exponent))]

minRegamma____minImgamma = ...
 [min(real(error_exponent)), min(imag(error_exponent))]

UREFcoarse = UREF(1:2^(kk-1):end);     % this uses the last used value of kk, i.e., kk=3
for kk = 1: 3
    max_rel_error(kk) = max( (abs( UREFcoarse - Uend(kk,:)) ) )/max( UREFcoarse );
end
max_rel_error = max_rel_error








