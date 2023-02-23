close all 
clear, clc

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

