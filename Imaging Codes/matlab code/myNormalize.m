function [normalize] =myNormalize(A)
[r,c] = size(A);
 max2=max(max(A));
 min2=min(min(A));
for i=1:r 
    for j=1:c 
        normalize(i,j)=round((A(i,j)-min2).*255./(max2-min2)); 
    end
end
end
% 
% function [normalize] =myNormalize(A)
%  max2=max(max(A));
%  min2=min(min(A));
% lg = log10(A+1);
% normalize=round((lg-min2)/(max2-min2))*255; 
%  
% end