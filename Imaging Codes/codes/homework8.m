% filters
%  img = imread('Lenna.png'); 
%  gray = rgb2gray(img);
% 
%  h = fir1(4,0.25,'low');
%  h = h'*h;
%  freqz2(h)
%  y=conv2(h,gray,'same')
% %  y = filter2(h,gray,'same');
%  figure, imshow(gray);
%  figure, imshow(uint8(y));


A = imread('fred.jpg');
A = imresize(A, [200 200]);
% A = im2double(A);
A = rgb2gray(A);
size(A)
subplot(2,2,1)
imagesc(A); title('original')
colormap('Gray')

M1 = size(A,1); N1 = size(A,2);M2 = size(A,1); N2 = size(A,2); 
DFA = fft2(double(A), M1, N1);

%......low-pass filtering by DFT windows......%
W1 = 100; W2 = 100;
r1 = zeros(M1,1);
r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
r2 = zeros(N2, 1);
r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
size(r1),size(r2)
w = r1*r2';
C = real ( ifft2 (DFA.*w));
figure(1);subplot(2,2,2)
imagesc(C); colormap('Gray')
title('W1 = W2 =40')


% W1 = 30; W2 = 30;
% r1 = zeros(M1,1);
% r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
% r2 = zeros(N2, 1);
% r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
% w = r1*r2';
% C = real ( ifft2 (w.*DFA));
% subplot(2,2,3)
% imagesc(C); colormap('Gray')
% title('W1= W2 =30')
% 
% W1 = 20; W2 = 20;
% r1 = zeros(M1,1);
% r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
% r2 = zeros(N2, 1);
% r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
% w = r1*r2';
% C = real ( ifft2 (w.*DFA));
% subplot(2,2,4)
% imagesc(C); colormap('Gray')
% title('W1= W2 =20')
% 
% %...low-pass filtering in spatial domain......%
% W=1;
% L = zeros(3);
% L(2,2) = 1/(2*W+1)^2;
% % L = [-1 -2 -1; -2 12 -2; -1 -2 -1]/16;
% filteredImage = imfilter(single(A), L);
% figure(2);subplot(231);imshow(filteredImage, [])
% title('W=1', 'FontSize', 12);
% 
% W=3;
% L = zeros(3);
% L(2,2) = 1/(2*W+1)^2;
% filteredImage = imfilter(single(A), L);
% subplot(232);imshow(filteredImage, [])
% title('W=3', 'FontSize', 12);
% 
% 
% W=8;
% L = zeros(3);
% L(2,2) = 1/(2*W+1)^2;
% filteredImage = imfilter(single(A), L);
% subplot(233);imshow(filteredImage, [])
% title('W=8', 'FontSize', 12);
% 
% 
% %......high-pass filtering by DFT windows......%
% W1 = 100; W2 = 100;
% r1 = zeros(M1,1);
% r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
% r2 = zeros(N2, 1);
% r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
% w = 1 - r1*r2';
% C = real ( ifft2 (DFA.*w));
% figure(3); subplot(231)
% imagesc(C); colormap('Gray')
% title('W1 = W2 =40')
% subplot(234); imshow(fftshift(w))
% title('h (normalized)')
% 
% 
% W1 = 30; W2 = 30;
% r1 = zeros(M1,1);
% r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
% r2 = zeros(N2, 1);
% r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
% w = 1 - r1*r2';
% C = real ( ifft2 (w.*DFA));
% subplot(232)
% imagesc(C); colormap('Gray')
% title('W1= W2 =30')
% subplot(235); imshow(fftshift(w))
% title('h (normalized)')
% 
% 
% W1 = 20; W2 = 20;
% r1 = zeros(M1,1);
% r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
% r2 = zeros(N2, 1);
% r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
% w = 1 - r1*r2';
% C = real ( ifft2 (w.*DFA));
% subplot(233)
% imagesc(C); colormap('Gray')
% title('W1= W2 =20')
% subplot(236); imshow(fftshift(w))
% title('h (normalized)')
% 
% %...high-pass filtering in spatial domain......%
% 
% W=1;
% L = -1 * zeros(3);
% L(2,2) = 1/2*(W+1)^2;
% % L = [-1 -2 -1; -2 12 -2; -1 -2 -1]/16;
% filteredImage = imfilter(single(A), L);
% figure(4);subplot(231);imshow(filteredImage, [])
% title('W=1', 'FontSize', 12);
% 
% W=3;
% L = -1 * zeros(3);
% L(2,2) = 1/2*(W+1)^2;
% filteredImage = imfilter(single(A), L);
% subplot(232);imshow(filteredImage, [])
% title('W=3', 'FontSize', 12);
% 
% W=8;
% L = -1 * zeros(3);
% L(2,2) = 1/2*(W+1)^2;
% filteredImage = imfilter(single(A), L);
% subplot(233);imshow(filteredImage, [])
% title('W=8', 'FontSize', 12);


%......band-pass filtering by DFT windows......%


%...band-pass filtering in spatial domain......%

