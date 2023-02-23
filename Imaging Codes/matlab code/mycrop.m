function [B] =mycrop(A, n1, n2, N1, N2)
% mycrop.m
% [B] =mycrop(A; n1; n2; N1; N2)
% crops the image A from location n1; n2
% with size N1; N2
for k = 0 : N1-1
    for l = 0 : N2-1
        B(k + 1, l + 1) = A(n1 + k + 1, n2 + l + 1);
    end
end