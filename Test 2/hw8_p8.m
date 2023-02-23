% Homework 8 Problem 8
% Solving BVP  
%  y'' = (y^2)/(2+x),  y(0)=1,  y(2)=1
%  by the Modified Picard method.

clear all

Solution_type = 'shallow';       % use "shallow" for the shallow solution and
                            %     "deep" for the deep solution, and
                            %     "deeper" for the deeper solution

% Set up the x-vector and the boundary conditions:
h = 0.1;
hsq = h^2;
x = 0 :h: 2;
xbvp = x(2:end-1);
N = length(xbvp);
alpha = 1;  % y(0)
beta = 1;   % y(2)
tolerance = 10^(-6);
max_Niterations = 1000;        % max allowed number of iterations


Cparam = input('enter the value for parameter C:  C = ');

% Build the diaginal and the 2 subdiagonals of matrix A:
a = 1*ones(1,N-1);
c = 1*ones(1,N-1);
b = -2*ones(1,N) - hsq*Cparam;

% Initial guess:
% Initial guess:
if strcmp(Solution_type, 'shallow') == 1
    Y = (alpha + beta)/2 * ones(1,N);
elseif strcmp(Solution_type, 'deep') == 1
    Nhalf = round(N/2);
    Y(1 : Nhalf) = 1 - 16*xbvp(1 : Nhalf);
    Y(Nhalf+1 : N) = -31 + 16*xbvp(Nhalf+1 :end);
elseif strcmp(Solution_type, 'deeper') == 1
    Nhalf = round(N/2);
    Y(1 : Nhalf) = 1 - 20*xbvp(1 : Nhalf);
    Y(Nhalf+1 : N) = -39 + 20*xbvp(Nhalf+1 :end);  
end

% Assign the initial error:
D(1)=1;
% # of iterations after which the solution is plotted within the loop below
Nplot = 20;

k=1;
while D(k) > tolerance && k < max_Niterations
    
%     Ysaved(k,:) = Y; % This line is needed only if I want to save iterations to look at them later.
    r = hsq * Y.^2./(2+xbvp) - hsq*Cparam*Y;
    r(1) = r(1) - alpha;
    r(N) = r(N) - beta;
    
    Ynew = thomas(a,b,c,r);
    D(k+1) = norm( Ynew - Y,inf);
    Y = Ynew;
    
    % Plot solutions after Nplot iterations:
    if rem(k, Nplot) == 0
        figure(80902);
        plot(x, [alpha Y beta]);
        xlabel('X','fontsize',16);  
        ylabel('Y','fontsize',16); 
        title(['iteration # ' int2str(k)], 'fontsize', 16)
        disp('Hit any key to continue with iterations.')
        pause
    end

    k=k+1;    
end

if k > max_Niterations
    disp(' !!!!!!!!!!!!!!!!!      Iterations do NOT converge         !!!!!!!!!!!!!')
    disp('  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ')
else
    sprintf('converges in %d  iterations', k)
    sprintf('last error = %d', D(k))
end

figure(80901);
subplot(2,1,1);
plot(x, [alpha Y beta]);
xlabel('X','fontsize',16);  
ylabel('Y','fontsize',16); 
subplot(2,1,2);
plot([1:k], log10(D),'r--');
xlabel('iteration`s number','fontsize',16);  
ylabel('log_{10} of error','fontsize',16); 
