function [y, I]=maxk(A,k)
[A,I]=sort(A,'descend');
y=A(1:k);
I=I(1:k);
end