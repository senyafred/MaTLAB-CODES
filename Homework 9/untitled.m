%This is the code for problem 1 hw 9
clear all
y0=0;
yf=1;
M = input('enter M = ');
x0=0;
xf=1;
h=1/M;
A= [M-1; M-1];


for k=1:M-1

   x(k)= x0 + h*(k);
   R(k)= ( 2*x(k) - 4 )/( 1 + x(k) )^2;

    for j=1:M-1
        A(k,j) = -((j*pi)^2) * sin(j*pi*x(k)) + ( -2/( 1 + x(k))^2 ) * sin(j*pi*x(k));
    end 

 
end
   
c= A\transpose(R);

for k=1:M-1

    Y(k)=0;
    for j=1:M-1
        Y(k)= Y(k)+c(j)*sin(j*pi*x(k));
    end
    Y(k)=Y(k)+x(k);

end

Y=[y0,Y,yf];

x=[x0,x,xf];

for i=1:length(x)
   y(i)= 2 * x(i)/( 1 + x(i));
end    

figure(2)
plot(x,Y,'r--')
hold on
plot (x,y)
hold off

figure(1)
plot(x,Y-y)