% LU decomposition
% function [L, U] = Gaussian_Elimination(A, b) 
% In fact this can be used to find an LU decomposition of A
function [L, U, x] = LU(A, b)
% ==== Description ====
% Solves a linear system Ax=b
% using Gaussian Elmination
%
% ==== input ====
% A ...... n-by-n matrix
% b ...... n-vector
%
% ==== output ====
% x ...... solution to Ax=b

verbose = false; % turn of this option to suppress outputs
% create an augmented matrix
M = [A, b];
if verbose
    disp('Initially, ');
    disp(M);
end

[n, ~] = size(A); % 
L = eye(n); % for LU decomposition
for i=1:n % for each row
    for j= i+1:n % for every other row below the current row
        m = M(j,i)/M(i,i); % multiplier
        L(j,i) = m; % for LU decomposition
        M(j,:) = M(j,:) - m * M(i,:);
        if verbose
            disp([ 'i = ', num2str(i), ', j = ', num2str(j)]);
            disp(M);
        end
    end
end
U = M(:,1:n); % for LU decomposition

% Now the left n-by-n part of M is upper triangular
x = BackSubstitution(M(:,1:n), M(:,end));
if verbose
    disp('solution = ');
    disp(x);
    disp('check; A*x =  ');
    disp(A*x);
    
end


end