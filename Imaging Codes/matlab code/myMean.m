function [Totalmean] =myMean(A)
[r,c] = size(A);
Totalmean = sum(A(:)/(r*c));
%meanIntensity = mean(img(:));
end