% % T = maketform('affine',[.5 0 0; .5 2 0; 0 0 1]);
% % tformfwd([10 20],T);
% % I = imread('cameraman.tif');
% % transformedI = imtransform(I,T);
% % % figure, imshow(I), figure, imshow(transformedI)
% 
% % src: http://www.imm.dtu.dk/~jmca/02501/exercises_old/image_warping/html/
% % load input image
% I = double(imread('fred.jpg'));
% [h,w,d] = size(I);
% 
% % show input image
% figure; image(I/255); axis image;
% title('input image');
% 
% % make transformation matrix (T)
% s  = 1;
% a  = pi/6;
% tx = 0;
% ty = 0;
% 
% T = [1 0 0; 0 0 0; -50 0 1];
% % T  = [ s*cos(a) s*sin(a) tx ; -s*sin(a) s*cos(a) ty ; 0 0 1 ];
% % T  = [ s*cos(a) s*sin(a) tx ; -s*sin(a) s*cos(a) ty ; 0 0 1 ];
% 
% % warp incoming corners to determine the size of the output image
% cp = T*[ 1 1 w w ; 1 h 1 h ; 1 1 1 1 ]; 
% Xpr = min( cp(1,:) ) : max( cp(1,:) ); % min x : max x
% Ypr = min( cp(2,:) ) : max( cp(2,:) ); % min y : max y
% [Xp,Yp] = ndgrid(Xpr,Ypr);
% [wp,hp] = size(Xp); % = size(Yp)
% 
% % do backwards transform (from out to in)
% n = wp*hp;
% X = T \ [ Xp(:) Yp(:) ones(n,1) ]';  % warp
% 
% % re-sample pixel values with bilinear interpolation
% clear Ip;
% xI = reshape( X(1,:),wp,hp)';
% yI = reshape( X(2,:),wp,hp)';
% Ip(:,:,1) = interp2(I(:,:,1), xI, yI, '*bilinear'); % red
% Ip(:,:,2) = interp2(I(:,:,2), xI, yI, '*bilinear'); % green
% Ip(:,:,3) = interp2(I(:,:,3), xI, yI, '*bilinear'); % blue
% 
% % show the warping result
% figure; image(Ip/255); axis image;
% title('warped image');
% 
% 
% % clear; 
% % close all;
% % img = imread('anas.jpg');
% % m =[1, -0.1, 0
% %     0, 1, 0.2
% %     0, 0, 1];
% % sz = size(img);
% % [X,Y]= meshgrid([1:sz(2)], [1:sz(1)]);
% % A = [reshape(Y,1,[]);
% %      reshape(X,1,[]);
% %      ones(1,length(Y(:)));];
% % AA = m*A;
% % AA = AA./[AA(3,:);AA(3,:);AA(3,:)];
% % AA = int32(AA);
% % a1 = AA(1,:) - min(AA(1,:)) +1;
% % a2 = AA(2,:) - min(AA(2,:)) +1;
% % 
% % 
% % szNew = sz;
% % szNew(1) = max(a1);
% % szNew(2) = max(a2);
% % 
% % ind = a1 + (a2-1)*szNew(1);
% % indOld = A(1,:) + (A(2,:)-1)*sz(1);
% % 
% % imgNew = uint8(zeros(szNew(1),szNew(2),3))*255;
% % 
% % imgNew(ind) = img(indOld);
% % imgNew(ind + szNew(1)*szNew(2)) = img(indOld + sz(1)*sz(2));
% % imgNew(ind + szNew(1)*szNew(2)*2) = img(indOld + sz(1)*sz(2)*2);
% % 
% % figure('position',[200,500,500,300]);
% % imshow(img);
% % figure('position',[800,500,500,300]);
% % imshow(imgNew);
% 
% 
% %%
% % 1. Implement the pseudo random noise dithered quantizer with
% % ∆ = 32, 64, 128. In each case try to choose the “best” possible
% % noise parameter A.
% 
% % A = imread('Lenna.png');
% % A = rgb2gray(A);
% % 
% % delta1 = 128;    
% % a = 512;                    % noise parapmeter
% % A1 = uint8(round(512*rand(a)));
% % A = A1+A;
% % figure; imshow(A,[])
% % B = delta *floor(A/delta) + delta/2;  % perform uniform quantization
% % figure; imshow(B,[])
% % % A = imnoise(A,'salt & pepper', 0.02);
% % 
% % delta1 = 32;    
% % delta2 = 64;
% % delta3 = 128;
% % a = 0.02;       % noise parameter



%%
% 2. Implement the warping functions for translation, rotation, wave,
% warp, swirl and glass. Try to use different parameters. My
% parameters x 0 , y 0 , etc., are set for a 512 × 512 image.


img = imread('fred.jpg');
A = rgb2gray(img);

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







%% 
% 3. Implement the median filter for your image. First randomly pick
% 4000 pixels and set their values to zero to obtain a pepper noise
% corrupted image. Then apply the median filter with W = 3.
% Repeat with 40000 pixels corrupted. What can you do to
% improve results now? (Hint: Try increasing W and/or running
% several iterations of the median filter.)
img = imread('Lenna.png');
img = rgb2gray(img);

[M,N] = size(img);

%% prepare a pepper noise corrupted image
% generate random binary matrix

f = 4000/(M*N);
noise = ones(M,N); % pre-allocate result
k = round(f*M*N); % number of 1's to place in result
noise(randperm(M*N,k)) = 0;

% zero out portion which coincides with zeros
A = double(img).*noise;
% B = medfilt2(A);
Z = zeros(M,N);
% set W
W = 1;

B = zeros(M,N);

for i=1:M
    for j=1:N
        x = (i-W):(i+W);
        y = (j-W):(j+W);
        B(i,j) = median(median(A(x(x>0 & x<=M),y(y>0 & y<=N))));
    end
end
imshow(uint8(B));



