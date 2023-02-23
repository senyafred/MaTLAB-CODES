function [Gx] =mySimple(bw)
 
bw=uint8((1/3)*(double(bw(:,:,1))+double(bw(:,:,2))+double(bw(:,:,3))));
bw=im2bw(bw);
[r,c]=size(bw);
bwdbl= double(bw);
%%%%Apply a mask 
maskx=[-1 0 1];
OUT = zeros(r-3,c-3);
for idx = 1:(r-3)
    for jdx=1:(c-3)
        bwsquare = bwdbl(idx:(idx+2),jdx:(jdx+2));
        res=maskx.*bwsquare;
        OUT(idx,jdx)=sum(sum(res));
    end
end
Gx = OUT;
end