function [B] =myConvolve(A, k)
% [r,c] = size(A);
% [m,n] = size(k);
% h = rot90(k, 2);
% center = floor((size(h)+1)/2);
% 
% Rep = zeros(r + m*2-2, c + n*2-2);
% for x = m : m+r-1
%     for y = n : n+r-1
%         Rep(x,y) = A(x-m+1, y-n+1);
%     end
% end
% B = zeros(r+m-1,n+c-1);
% for x = 1 : r+m-1
%     for y = 1 : n+c-1
%         for i = 1 : m
%             for j = 1 : n
%                 B(x, y) = B(x, y) + (Rep(x+i-1, y+j-1) * h(i, j));
%             end
%         end
%     end
% end

[r,c] = size(A);
[m,n] = size(k);
E = rot90(k,2);
center = floor((size(E)+1)/2);
left = center(2)-1;
right = n - center(2);
top = center(1)-1;
bottom = m - center(1);
Rep = zeros(r+top+bottom,c+left+right);
for x= 1+top : r+top
    for y= 1+left : c+left
        Rep(x,y)= A(x-top,y-left);
    end
end
WW = zeros(r,c);
for x = 1:r
    for y = 1:c
        for i = 1:m
            for j = 1:n
                q = x-1;
                w = y-1;
                WW(x,y)= WW(x,y)+ (Rep(i+q,j+w))* E(i,j);
           end
       end
   end
end