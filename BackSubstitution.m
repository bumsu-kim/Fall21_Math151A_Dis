function x = BackSubstitution(U, b)
% ==== Description ====
% Solves a upper triangluar system Ux=b
%
% ==== input ====
% U ...... n-by-n UPPER TRIANGULAR matrix
% b ...... n-vector
%
% ==== output ====
% x ...... solution to Ux=b

n = length(b);
x = zeros(n,1);
for j = n:-1:1 % NOTE - n:1 won't work
    x(j) = (b(j) - dot(U(j, j+1:end), x(j+1:end)))/U(j,j);
end

end