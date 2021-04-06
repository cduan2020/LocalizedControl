function [y, I]=mink(A,k)
% n=size(A,1);
% A=A-eye(n);
[A,I]=sort(A,1);
y=A(1:k,:);
I=I(1:k,:);
end