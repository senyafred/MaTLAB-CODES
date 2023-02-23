% Fourier_Transform.pdf Assignment

%% 1. Calculate the DFT of your image. Show the magnitude and phase
%     both before and after using fftshift. Use the log10 point function
%     on magnitude plots and normalize as necessary.

% A = imread('anas.jpg');
% A = A(:,:,1);
% % figure
% imshow(A);
% 
% % DFT and Plot the magnitude and phase in degrees, before fftshift:
% y1 = fft2(A);                               % Compute DFT of x
% m1 = abs(A);                                % Magnitude
% p1 = unwrap(angle(y1));                      % phase
% f = (0:length(y1)-1)*100/length(y1);         % Frequency vector
% 
% subplot(2,1,1)
% plot(f,m1)
% title('Magnitude')
% ax = gca;
% ax.XTick = [15 40 60 85];
% 
% subplot(2,1,2)
% plot(f,p1*180/pi)
% title('Phase')
% ax = gca;
% ax.XTick = [15 40 60 85];
% 
% % DFT and Plot the magnitude and phase in degrees,  fftshift:
% y2 = fftshift(y1);                         % Compute DFT of x
% m2 = abs(A);                                % Magnitude
% p2 = unwrap(angle(y2));                     % phase
% f = (0:length(y2)-1)*100/length(y2);        % Frequency vector
% 
% figure
% subplot(2,1,1)
% plot(f,m2)
% title('Magnitude')
% ax = gca;
% ax.XTick = [15 40 60 85];
% 
% subplot(2,1,2)
% plot(f,p2*180/pi)
% title('Phase')
% ax = gca;
% ax.XTick = [15 40 60 85];
% 
% % log10 point function on magnitude plots and normalize as necessary
% m11 = log10(abs(double(A)));
% figure
% plot(f,m2)
% title('log10 point function on magnitude plot')
% ax = gca;
% ax.XTick = [15 40 60 85];

%% 2)
% important: https://stackoverflow.com/questions/28818632/subsample-an-image-using-a-for-loop
% 2. Subsample your image by 2 in each direction and calculate the
% DFT of the result in two ways:
% 2.1 Since the subsampled image has half the dimensions of the
% original, calculate the DFT to the point of the reduced dimensions.
% 2.2 Calculate the DFT to the point of the original dimensions.
% Show magnitude and phase plots for both. Compare the results
% to 1. above. Explain the differences?

% sub sampling image by a factor of 2 
I=imread('anas.jpg');
I = I(:,:,1);
J = imresize(I, 0.5);
% J=downsample(I,2);
% J=downsample(J',2)';
whos I J
imshow(J)
figure, imshow(I)

% 2.1: calculating the DFT to the point of the reduced dimensions.
DFT_J = fft2(J);                               % Compute DFT of x
MAG_J = abs(J);                                % Magnitude
PHASE_J = unwrap(angle(DFT_J));               % phase
f1 = (0:length(DFT_J)-1)*100/length(DFT_J);    % Frequency vector

subplot(2,1,1)
plot(f1,MAG(J))
title('Magnitude')
ax = gca;
ax.XTick = [15 40 60 85];

subplot(2,1,2)
plot(f1,PHASE_J*180/PHASE_J)
title('Phase')
ax = gca;
ax.XTick = [15 40 60 85];


% 2.2: Calculate the DFT to the point of the original dimensions.
DFT_I = fft2(I);                               % Compute DFT of x
MAG_I = abs(I);                                % Magnitude
PHASE_I = unwrap(angle(DFT_I));               % phase
f2 = (0:length(DFT_I)-1)*100/length(DFT_I);    % Frequency vector

subplot(2,1,1)
plot(f2, MAG_I)
title('Magnitude')
ax = gca;
ax.XTick = [15 40 60 85];

subplot(2,1,2)
plot(f2,PHASE_I*180/PHASE_I)
title('Phase')
ax = gca;
ax.XTick = [15 40 60 85];

% smallImage = I(1:2:end, 1:2:end);  % Remember (row, column) = (y, x) NOT (x,y)
% imshow(smallImage)
% F = fft2(A);
% figure
% imagesc(abs(F));
% figure
% imagesc(abs(fftshift(F)))
% figure
% imagesc(log10(abs(F)+1))
% figure
% imagesc(log10(abs(fftshift(F))+1))

%%
% 3. Calculate the convolution of your image with itself by using DFTs.
% M1 = size(A,1);
% N1 = size(A,2);
% DFA = fft2(A, M1+M1-1,N1+N1-1);
% DFB = fft2(A, M1+M1-1,N1+N1-1);
% DFC = DFA.*DFB;
% imagesc(log10(abs(DFC)))
% imagesc(log10(abs(fftshift(DFC))+1))
%%
% 4. Do the processing Pages 19-20. Show the resulting images and
% the MSE plot.
% A1 = imread('Lenna.png');
% A1 = rgb2gray(A1);
% DFA = fft2(A, M1, N1);
% DFA2 = fftshift(DFA);
% 
% figure
% subplot(2,3,1)
% imagesc(A1)
% colormap('Gray')
% 
% M1 = size(A1,1); N1 = size(A1,2);M2 = size(A1,1); N2 = size(A1,2); 
% W1 = 1; W2 = 1;
% r1 = zeros(M1,1);
% r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
% r2 = zeros(N2, 1);
% r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
% w = r1*r2';
% DFC = DFA.*w;
% C = real ( ifft2 (w.*DFA));
% subplot(2,3,2)
% imagesc(C); colormap('Gray')
% title('W1= W2 =1')
% 
% W1 = 10; W2 = 10;
% r1 = zeros(M1,1);
% r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
% r2 = zeros(N2, 1);
% r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
% w = r1*r2';
% C = real ( ifft2 (w.*DFA));
% subplot(2,3,3)
% imagesc(C); colormap('Gray')
% title('W1 = W2 =10')
% 
% 
% W1 = 20; W2 = 20;
% r1 = zeros(M1,1);
% r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
% r2 = zeros(N2, 1);
% r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
% w = r1*r2';
% C = real ( ifft2 (w.*DFA));
% subplot(2,3,4)
% imagesc(C); colormap('Gray')
% title('W1= W2 =20')
% 
% W1 = 30; W2 = 30;
% r1 = zeros(M1,1);
% r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
% r2 = zeros(N2, 1);
% r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
% w = r1*r2';
% C = real ( ifft2 (w.*DFA));
% subplot(2,3,5)
% imagesc(C); colormap('Gray')
% title('W1= W2 =30')
% 
% W1 = 40; W2 = 40;
% r1 = zeros(M1,1);
% r1(1:W1 + 1,M1 - W1 + 1 :M1) = 1;
% r2 = zeros(N2, 1);
% r2(1 : W2 + 1, N1 - W2 + 1 : N1) = 1;
% w = r1*r2';
% C = real ( ifft2 (w.*DFA));
% subplot(2,3,6)
% imagesc(C); colormap('Gray')
% title('W1= W2 =40')
% 
% % ploting mean square error of Lenna image
% size(A1)
% size(C)
% MSE = immse(DFA, DFC)