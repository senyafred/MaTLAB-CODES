J= imread('fred.jpg');
% %[M,N]= size(Pic);
% %imtool(Pic);
 J = rgb2gray(J);
% %imshow(I);
% [m,n]=size(I);
% training_set= double(I);
% training_set = reshape(training_set, m*n,1);
% s = 5:-1:-1;
%       len = 2.^s;
%       [partition, codebook] = lloyds(training_set, len);
%       [PicLloyd ,index] = imquantize(I,partition,codebook);
% 
% % figure();
% % x = uint8(PicLloyd);
% %  imshow(x);
% % err = immse(I,x);
% % fprintf('\the mean-squared error is %0.4f\n',err);
% %max_t=double(max(gray(:)));
% % PSNR calculation
% %squaredErrorImage = (double(gray) - double(PicLloyd)) .^ 2;
% %mse = sum(sum(squaredErrorImage)) / (M * N);
% %PSNR = 20 * log10(max_t/ mse)
% 
% 
%  %I = [0 4 4;
%   %     -2 1 3;
%   %      5 1 1];
%   tic
%  k = [-1 0 1;
%        -2 0 2;
%        -1 0 1];
% Hsl = abs(myConvolve(double(I), double(k)));
% %Hsl = conv2(I, I);
% display(Hsl);
% % display(Hsl1)
% Hsl=myNormalize(Hsl);
% figure();
% imshow(Hsl);
% toctic;

% imtool(J);
HH=[-1 0 1; -2 0 2 ; -1 0 1];
% II = (convolve(double(J),HH))));
II = abs(myConvolve(double(J),double(HH)));
sca11 = round(II-min(II(:)))*255 ./ (max(II(:)-min(II(:))));
RR = uint8(sca11);
% MM = imshow(RR)
toc;

% HH=[-1 0 1; -2 0 2 ; -1 0 1];
% B = imread('kofi.png');
% imtool(B);
% II = rgb2gray(B);
% RR1 = round(log10(abs(conv2(double(II),HH,'same'))));
RR1 = abs(conv2(double(HH),double(J)));
sca111 = (RR1-min(RR1(:)))*255 ./ (max(RR1(:)-min(RR1(:))));
RRR = uint8(sca111);
% MMM= imshow(RRR)
% MMM = uint8(RR);
% imagesc(MMM);
toc;
subplot(2,2,1),imshow(J);
title('REAL IMAGE');
subplot(2,2,2),imshow(RR);
title('CONV FUNCTION');
subplot(2,2,3),imshow(RRR);
title('CONV2');





