function [Totalvar] =myVariance(A)
[r,c] = size(A);
 Totaldiff=(A-myMean(A)).^2;
 Totalsum=sum(Totaldiff(:));
 nele=(r*c)-1;
 Totalvar=Totalsum/nele;
end