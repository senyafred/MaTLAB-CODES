% % I = imread('coins.png');
% % imshow(I)
% % imtool(I)
% P=imread('peppers.png');
% imtool(P)
% m=fspecial('motion',50,50);
% I=imfilter(P,m);
% imtool(I)

%  I=imread('peppers.png');
%  imshow(I)
% % size(I)
% % impixel(I,20,40)
% % [b,bmap]=imread('emu.tif');
% A=rgb2gray(I);
% imshow(A)
% 
% c=imread('cameraman.tif');
% od=double(c);
% imshow(od/255)
% a = c>120;  
% imshow(a)

% RESIZE AN IMAGE
% n=imread('cameraman.tif');
% [rows,columns]=size(n);
% c=zeros(rows/2,columns/2);
% j=1;i=1;
% for x=1:2:rows
%     for y=1:2:columns
%         c(i,j)=n(x,y);
%         j=j+1;
%     end
%     i=i+1;
%     j=1;
% end
% figure, imshow(n);
% figure, imshow(c/255);
% figure,imagesc(c),colormap(gray)
    

%  W = imread('cameraman.tif');
%  id1 = double(W);
% % id1=id+100;
% % imshow(uint8(id1))
% id2=id1-50;
% imshow(uint8(id2))

% u=imread('peppers.png');
% ad=im2double(u);
% x=ad;y=ad;
% [r,c]=size(ad);
% factor=1;
% for i=1:r
%     for j=1:c
%         x(i,j)=factor*log(1+ad(i,j));%Reduce brightness
%         y(i,j)=factor*ad(i,j)^2;%Enhances brightness
%     end
% end
% subplot(2,2,1);imshow(ad);title('Before')
% subplot(2,2,2);imshow(x);title('After')
% subplot(2,2,3);imshow(ad);title('Before')
% subplot(2,2,4);imshow(y);title('After')

% I =imread('pout.tif');
% figure;
% subplot(1,2,1);imshow(I);
% subplot(1,2,2);imhist(I);
% 
% imh=imadjust(I,[0.3,0.6],[0.0,1.0]);
% P=histeq(I);
% subplot(2,2,1);imshow(P);
% subplot(2,2,2);imhist(P);
% subplot(2,2,3);imshow(imh);
% subplot(2,2,4);imhist(imh);

% % Low pass filters
% R = imread('cameraman.tif');
% f=ones(3)/9;
% R1=filter2(f,R,'same');
% figure;
% imshow(R1/255);
% % OR
% R2=fspecial('average',[3,3]);
% I2=filter2(R2,R,'valid');
% figure;imshow(I2/255);

% W=imread('trees.tif');
% an=imnoise(W,'gaussian',0.01);
%  %remove Gaussian noise using Gaussian Filter
% sigma=3;
% cutoff=ceil(3*sigma);
% h=fspecial('gaussian',2*cutoff+1,sigma);
% out=conv2(an,h,'same');
% %Using wiieners
% W1=wiener2(an,[5 5]);
% out1=conv2(W,h,'same');
% subplot(2,2,1);imshow(an);
% subplot(2,2,2);imshow(out/255);
% subplot(2,2,3);imshow(W);
% subplot(2,2,4);imshow(out1/255);
% figure;imshow(W1);
% figure;surf(1:2*cutoff+1,1:2*cutoff+1,h);


% y=imread('cameraman.tif');
% isp=imnoise(y,'salt & pepper',0.1);
% figure;imshow(isp);
% %using average filter
% f=fspecial('average');
% I=filter2(f,isp);
% imshow(uint8(I));
% %using median filter
% med=medfilt2(isp);
% figure;imshow(med);



% 
% %Histogram equalization without function
% T=imread('cameraman.tif');
% [r,c]=size(T);
% ah=uint8(zeros(r,c));
% n=r*r;
% f=zeros(256,1);
% pdf=zeros(256,1);
% cdf=zeros(256,1);
% cum=zeros(256,1);
% out=zeros(256,1);
% for i=1:r
%     for j=1:c
%         value=T(i,j);
%         f(value+1)=f(value+1)+1;
%         pdf(value+1)=f(value+1)/n;
%     end
% end
% sum=0; L=255;
% for i=1:size(pdf)
%     sum=sum + f(i);
%     cum(i)=sum;
%     cdf(i)=cum(i)/n;
%     out(i)=round(cdf(i)*L);
% end
% for i=1:r
%     for j=1:c
%         ah(i,j)=out(T(i,j)+1);
%     end
% end
% figure,imshow(ah);
% he=histeq(T);
% figure,imshow(he);title('using function')


% %Implimenting Edge
% ic=imread('circuit.tif');
% px=[-1 0 1;-1 0 1;-1 0 1];
% icx=filter2(px,ic);
% figure,imshow(icx/255);
% py=px';
% icy=filter2(py,ic);
% figure,imshow(icy/255);
% pedge=sqrt(icx.^2 + icy.^2);
% figure,imshow(pedge/255);
% fe=im2bw(pedge/255,0.3);
% figure,imshow(fe);

 

% % EDGE DETECTION
% E = imread('coins.png');
% imshow(E);
% BW1=edge(E ,'canny');
% BW2=edge(E,'sobel');
% subplot(1,2,1);imshow(BW1);
% subplot(1,2,2);imshow(BW2);

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%Edge Detection algorithm
% original=imread('q1jjy.png');
% %imshow(original);
% clipped=original(:,1:550,:);
% %imshow(clipped);
% %%%%%%%convert the image to black and white
%  %bw=uint8((1/3)*(double(clipped(:,:,1))+double(clipped(:,:,2))+double(clipped(:,:,3))));
%  bw=im2bw(clipped);
%  figure()
% imshow(bw);
% [r,c]=size(bw);
% bwdbl= double(bw);
% %%%%Apply a mask in x
% maskx=[-1 -2 -1;0 0 0;1 2 1];
% OUT = zeros(r-3,c-3);
% for idx = 1:(r-3)
%     for jdx=1:(c-3)
%        % hold on
%         %rectangle('Position',[idx jdx 3 3]);
%         %pause
%         bwsquare = bwdbl(idx:(idx+2),jdx:(jdx+2));
%         res=maskx.*bwsquare;
%         OUT(idx,jdx)=sum(sum(res));
%     end
% end
% Gx = OUT;
% figure()
% imshow(Gx)
% 
% %%%%Apply another mask in y
% masky=[-1 0 1;-2 0 2;-1 0 1];
% for idx = 1:(r-3)
%     for jdx=1:(c-3)
%        % hold on
%         %rectangle('Position',[idx jdx 3 3]);
%         %pause
%         bwsquare = bwdbl(idx:(idx+2),jdx:(jdx+2));
%         res=masky.*bwsquare;
%         OUT(idx,jdx)=sum(sum(res));
%     end
% end
% Gy = OUT;
% figure()
% imshow(Gy)
% 
% %%Normalize the result of both masks
% G = sqrt(Gx.^2 + Gy.^2);
% 
% figure(),imshow(G);


        







