function [A_unique, count] = unique_tol(A, tol)
%   [A_unique, count] = unique_tol(A, tol)
%   Return all unique values and counts in A up to tol

A = A(~isnan(A));
A_int = round(A / tol);
A_unique = unique(A_int);
count = histc(A_int, A_unique);
A_unique = A_unique * tol;

end

