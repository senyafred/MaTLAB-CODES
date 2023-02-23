% %%%Reading the image
bw = imread('senya.jpg');
I = rgb2gray(bw);

%%%QUESTION 1
%%%%%% SIMPLE
Gx = mySimple(bw);
% figure();
% subplot(1,2,1),imshow(I),title('Original Image');
% subplot(1,2,2),imshow(Gx),title('Simple Edge');

%%%%Prewit
BW1=edge(I,'prewitt','vertical');
BW2=edge(I,'prewitt','horizontal');
BW3=edge(I,'prewitt','both');

%%%%Sobel
SW1=edge(I,'sobel','vertical');
SW2=edge(I,'sobel','horizontal');
SW3=edge(I,'sobel','both');
% figure();
% subplot(1,3,1),imshow(BW1),title('Prewitt Vertical');
% subplot(1,3,2),imshow(BW2),title('Prewitt Horizontal');
% subplot(1,3,3),imshow(BW3),title('Prewitt Combined');
% figure();
% subplot(1,3,1),imshow(SW1),title('Sobel Vertical');
% subplot(1,3,2),imshow(SW2),title('Sobel Horizontal');
% subplot(1,3,3),imshow(SW3),title('Sobel Combined');


%%%QUESTION 2
gaussian=edge(I,'log');
% figure();
% imshow(gaussian),title('Laplacian operator');


%%QUESTION 3
[u,v] = size(I);
newimage=double(I)+10*randn(u,v);
P=edge(newimage,'log');
figure();
imshow(P),title('Noise operator');


