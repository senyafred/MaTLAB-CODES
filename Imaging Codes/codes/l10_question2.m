close all 
clear, clc

img = imread('fred.jpg');

A = rgb2gray(img);
A = imresize(A, [200 200]);

[M,N] = size(A);

%% Question 2

% TRANSLATION
B1 = zeros(M,N);
for k=1:M
    for l=1:N
%         prepare rotation warping functions        
        x = k+50;
        y = l;
%         check for range of validity        
        if ((1<=x) && (x<=M)) && ((1<=y)&&(y<=M))
            B1(k,l) = A(x,y);
        else
            B1(k,l) = 0;
        end
    end
end

% ROTATION
B2 = zeros(M,N);
% calculate the center of the image
center = size(A)/2 + .5;
theta = pi/6;
for k=1:M
    for l=1:N
%         prepare rotation warping functions
        x = (k-center(1))*cos(theta) + (l-center(2))*sin(theta)+center(1);
        y = -(k-center(1))*sin(theta)+(l-center(2))*cos(theta)+center(2);
%         round x and y to prevent index error
        x = round(x);
        y = round(y);
%         check for range of validity
        if ((1<=x) && (x<=M)) && ((1<=y)&&(y<=M))
            B2(k,l) = A(x,y);
        else
            B2(k,l) = 0;
        end        
    end
end

% WAVE
B3 = zeros(M,N);
for k=1:M
    for l=1:N
%         prepare rotation warping functions
        x = k+20*sin(((2*pi)/128)*l);
        y = l;
%         round x to prevent index error
        x = round(x);
%         check for range of validity
        if ((1<=x) && (x<=M)) && ((1<=y)&&(y<=M))
            B3(k,l) = A(x,y);
        else
            B3(k,l) = 0;
        end        
    end
end

% WARP
B4 = zeros(M,N);
% calculate the center of the image
center = size(A)/2 + .5;
for k=1:M
    for l=1:N
%         prepare rotation warping functions
        x = (sign(k-center(1))/center(1))*(k-center(1))^2 + center(1);
        y = l;
%         round x to prevent index error
        x = round(x);
        y = round(y);
%         check for range of validity
        if ((1<=x) && (x<=M)) && ((1<=y)&&(y<=M))
            B4(k,l) = A(x,y);
        else
            B4(k,l) = 0;
        end
    end
end

% SWIRL
B5 = zeros(M,N);
% calculate the center of the image
center = size(A)/2 + .5;
for k=1:M
    for l=1:N
%         calculate r which determines theta
        r = ((k-center(1))^2 + (l-center(2))^2)^.5;
        theta = (pi/512)*r;
%         prepare rotation warping functions
        x = (k-center(1))*cos(theta)+(l-center(2))*sin(theta) + center(1);
        y = -(k-center(1))*sin(theta)+(l-center(2))*cos(theta) + center(2);
%         round x to prevent index error
        x = round(real(x));
        y = round(real(y));
%         check for range of validity
        if ((1<=x) && (x<=M)) && ((1<=y)&&(y<=M))
            B5(k,l) = A(x,y);
        else
            B5(k,l) = 0;
        end
    end
end

% GLASS
B6 = zeros(M,N);
for k=1:M
    for l=1:N
%         prepare rotation warping functions
        x = k + (rand(1,1) - .5)*10;
        y = l + (rand(1,1) - .5)*10;
%         round x and y to prevent index error
        x = round(x);
        y = round(y);
%         check for range of validity
        if ((1<=x) && (x<=M)) && ((1<=y)&&(y<=M))
            B6(k,l) = A(x,y);
        else
            B6(k,l) = 0;
        end
    end
end