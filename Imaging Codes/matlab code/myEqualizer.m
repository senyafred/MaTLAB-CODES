function [ah] =myEqualizer(A)
[r,c]=size(A);
ah=uint8(zeros(r,c));
n=r*r;
f=zeros(256,1);
pdf=zeros(256,1);
cdf=zeros(256,1);
cum=zeros(256,1);
out=zeros(256,1);
for i=1:r
    for j=1:c
        value=A(i,j);
        f(value+1)=f(value+1)+1;
        pdf(value+1)=f(value+1)/n;
    end
end
sum=0; L=255;
for i=1:size(pdf)
    sum=sum + f(i);
    cum(i)=sum;
    cdf(i)=cum(i)/n;
    out(i)=round(cdf(i)*L);
end
for i=1:r
    for j=1:c
        ah(i,j)=out(A(i,j)+1);
    end
end
end





