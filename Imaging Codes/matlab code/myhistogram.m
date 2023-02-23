function [h] =myhistogram(A)
  h=zeros(256,1);
 for l=0:255
     h(l+1) = sum(sum(A==l));
 end
 bar(0:255,h);
end