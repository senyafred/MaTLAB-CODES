%%%Reading the image
I = imread('fred.jpg');
I = rgb2gray(I);
%figure();
%imshow(I),title('Original Image');

%%%Flipping the image
h=I(:,(end:-1:1),:);     %%Vertical
V=I((end:-1:1),:,:);     %%Horizontal
%D=double(I);
D=I';                    %%Diagonal
%figure();
%subplot(2,2,1),imshow(h),title('Vertical flip');
%subplot(2,2,2),imshow(V),title('horizontal flip');
%subplot(2,2,3),imshow(uint8(D)),title('Diagonal flip');


%%%converting the image
img=double(I);


%%%%Cropping the image
 B=mycrop(img,440,680,660,660);
 display(size(B));
 %figure();
 %imshow(uint8(B));

% %%Calculating mean
M=myMean(img);
display(M);
 
 %%Calculating the Variance
V = myVariance(img);
 display(V);
 
 %%The histogram of the image
% myhistogram(img);
 %title('Histogram of image');
 %display(size(img));
 
 %%Finding the square-root
 sqrtroot=img.^(.5);
 sqrtround = round(sqrtroot);
 %myhistogram(sqrtround);
 %title('Histogram of Square root image');

 
%%%Normalize
%normalize = myNormalize(sqrtround);
%myhistogram(normalize);
%title('Histogram of Square root/Normalize image');

%%%%%%QUESTION 5
T=myNormalize(sqrtroot);
%myhistogram(T);
%title('Histogram of Normalize image')

%%%%%QUESTION 6
[r,c] = size(img);
for i=1:r
    for j=1:c
        B(i,j) = log10(img(i, j) + 1);
    end
end
imgB=myNormalize(B);
%subplot(1,2,1),imshow(uint8(img));
%subplot(1,2,2),imshow(uint8(imgB));
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%QUESTION 1
% %The histogram of hb(L)=l
% 
% 
% %%%QUESTION 2
% R1 = [0, 30]; R2 = [30, 120]; R3 = [120, 150]; R4 = [150, 240];
% R5 =[240, 255];
% 
% B1 = 255 * ((img >= 0)&(img <= 30));
% B2 = 255 * ((img >= 30)&(img <= 120));
% B3 = 255 * ((img >= 120)&(img <= 150));
% B4 = 255 * ((img >= 150)&(img <= 240));
% B5 = 255 * ((img >= 240)&(img <= 255));
% 
% %imshow(uint8(B1));
% m1 = sum(sum(B1.*img))/sum(sum(B1));
% m2 = sum(sum(B2.*img))/sum(sum(B2));
% m3 = sum(sum(B3.*img))/sum(sum(B3));
% m4 = sum(sum(B4.*img))/sum(sum(B4));
% m5 = sum(sum(B5.*img))/sum(sum(B5));
% 
% C = m1*B1 + m2*B2 + m3*B3 + m4*B4 + m5*B5;
% C=myNormalize(C);
% %imshow(uint8(C));
% 
% %%%QUESTION 4
% %ah = myEqualizer(img);  my function
% ah=histeq(uint8(img));
% %subplot(2,2,1),imshow(uint8(img));
% %subplot(2,2,2),myhistogram(img);
% %subplot(2,2,3),imshow(uint8(ah));
% %subplot(2,2,4),myhistogram(ah);
% 
% %%%QUESTION 5
C = log10(img+1);
C=myNormalize(C);
original = histeq(I);
logimage = histeq(uint8(C));
% figure();
% subplot(2,2,1),imshow(I),title('Image');
% subplot(2,2,2),myhistogram(original),title('Histogram Equalization of Image');
% subplot(2,2,3),imshow(logimage),title('Processed Image');
% subplot(2,2,4),myhistogram(logimage),title('Processed Histogram');


%% Question 6
A = double(A);
   B = log10(A+1);
   mxb = max(max(B));mnb = min(min(B));
   B = round(((B-mn)/(mxb-mnb))*255);
   B = uint8(B);
   figure()
   subplot(2,2,1)
   imshow(uint8(A))
   title('Original image A')
   subplot(2,2,2)
   imhist(A)
   title('H_A(l)')
   xlabel('Intensity Values Pixels')
   ylabel('Frequency')
   subplot(2,2,3)
   imshow(B)
   title('B=log10(A(i,j)+1)')
   subplot(2,2,4)
   imhist(B)
   title('H_B(l)')
   xlabel('Intesity value Pixels')
   ylabel('Frequency')






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %f = @(sigma,x,MU)(1\sqrt(2*pi*SIGMA^2)*exp((-(x-MU)^2)/(2*SIGMA^2)));
% %T = @(x,f)(x*f);
% 
% %%%%%%%QUESTION 1
Mu = 2; 
sigma = 3;
A = normrnd(Mu, sigma, 256, 256);

%%%Finding the Mean
MeanA=myMean(A);
display(MeanA);

%%Finding the Variance
VarA=myVariance(A);
display(VarA);

%%Normalizing A
B=myNormalize(A);
% 
% %%%%Sample mean of B
 [N,M]= size(B);
 H = myhistogram(B);
% % Calculating mean of B:
 pmf = H./(N*M);
% MeanB=0;
% for l=0:255
%     MeanB=MeanB+(l*pmf(l+1)); 
% end
%  %%Calculating Variance of B:
% VarB=0;
% for l=0:255
%     VarB=VarB+((l-MeanB)^2*pmf(l+1)); 
% end
% %display(MeanB),display(VarB);
% figure(),imhist(B);
% title('Histogram of Normalized image');


%%%Histogram Matching
im1 = img;
im2 = B;
M = zeros(256,1,'uint8'); 
hist1 = myhistogram(im1); 
hist2 = myhistogram(im2);
cdf1 = cumsum(hist1) / numel(im1); %compute the cdf of im1
cdf2 = cumsum(hist2) / numel(im2); %compute the cdf of im2
 
for idx = 1 : 256
    diff = abs(cdf1(idx) - cdf2);
    [~,ind] = min(diff);
    M(idx) = ind-1;
end

%matching here
out = M(double(im1)+1);
% figure();
% subplot(2,3,1),imshow(uint8(im1));
% title('Original image');
% subplot(2,3,2),imshow(uint8(im2)),title('image B');
% subplot(2,3,3),imshow(uint8(out));
% title('matched image');
% subplot(2,3,4),imhist(uint8(im1));
% title('Histogram of original image');
% subplot(2,3,5),imhist(uint8(im2));
% title('Histogram of B');
% subplot(2,3,6),imhist(uint8(out));
% title('Histogram of matched image');

%%%%The point function
intensity = [0:255];
% figure();
% plot(intensity,M);
% xlabel('Intensity Values');
% ylabel('Differences');
% title('Matching Point Function');

%%%% UNDOING
A1 = log10(img+1);
A1=myNormalize(A1);
logimage = histeq(uint8(A1));
matchedimage = histeq(out);
%% subplot(2,3,1),imshow(logimage),title('Log Image');
%%subplot(2,3,4),myhistogram(logimage),title('Log Histogram');
%%subplot(2,3,2),imshow(matchedimage),title('matched Image');
%%subplot(2,3,5),myhistogram(matchedimage),title('Histogram of matched Image');
%%subplot(2,3,3),imshow(uint8(im2)),title('Image B');
%% subplot(2,3,6),myhistogram(uint8(im2)),title('Histogram of B');

% subplot(2,2,1),imshow(logimage),title('Log Image');
% subplot(2,2,2),myhistogram(logimage),title('Log Histogram');
% subplot(2,2,3),imshow(matchedimage),title('matched Image');
% subplot(2,2,4),myhistogram(matchedimage),title('Histogram of matched Image');


%%%%%QUESTION 4 QUANTIZATION
delta1=4;delta2=8;delta3=16;delta4=64;

Q1 = delta1*floor(img/delta1) + delta1/2;
Q2 = delta2*floor(img/delta2) + delta2/2;
Q3 = delta3*floor(img/delta3) + delta3/2;
Q4 = delta4*floor(img/delta4) + delta4/2;

%View image
% figure();
% subplot(2,2,1),imshow(uint8(Q1));
% subplot(2,2,2),myhistogram(uint8(Q1));
% subplot(2,2,3),imshow(uint8(Q2));
% subplot(2,2,4),myhistogram(uint8(Q2));
% figure();
% subplot(2,2,1),imshow(uint8(Q3));
% subplot(2,2,2),myhistogram(uint8(Q3));
% subplot(2,2,3),imshow(uint8(Q4));
% subplot(2,2,4),myhistogram(uint8(Q4));


%%Calculating MSQE
[w,x] = size(img);
MSQR1 = sum(sum(((img - Q1).^2)/(w*x)));
MSQR2 = sum(sum(((img - Q2).^2)/(w*x)));
MSQR3 = sum(sum(((img - Q3).^2)/(w*x)));
MSQR4 = sum(sum(((img - Q4).^2)/(w*x)));
% display(MSQR1);display(MSQR2);display(MSQR3);display(MSQR4);


%%%Companding %QUESTION 5
% [N,M]=size(img);
%  
% H = myhistogram(img);
% pmf = H./(N*M);
% 
% %pmfA=(histfreq(img)/n1);
% g1=zeros(256,1);
% for l=0:255
%     t=0;
%     for k=0:l
%         t=t+sum(pmf(k+1));
%     end
%   g1(l+1)=t;
% end
% 
% G = zeros(256/delta1,1,'uint8');
% 
% for k = 0: 255
%     for idx = 0 : (256/delta1)-1
%         diff = abs(255*g1(idx) - (((idx+1)*delta1) + (delta1/2)));
%         [~,ind] = min(diff);
%         G(idx) = ind-1;
%     end
% end
% display(G);
pA(l) = [.1, 0, .3, .2, 0, 0, .3, .1];
pB(k) = [.2, 0, 0, .1, .4, .3];



