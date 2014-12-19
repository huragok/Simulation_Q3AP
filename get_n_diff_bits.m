function B = get_n_diff_bits(Nbps)
%   B = get_n_diff_bits(Nbps)
%   Get the number of different bits between integers in [0, 2 ^ Nbps - 1]
% _____________________________________________________________________________
%	Inputs:
% 		Nbps:       scalar, number of bits per symbol
%	Outputs:
%		B:			2 ^ Nbps-by-2 ^ Nbps matrix, B(i, j) is the difference
%                   in bits between integer (i - 1) and (j - 1)
% _____________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 12/18/2014
% Codename: Dunkirk
% _____________________________________________________________________________

Q = 2 ^ Nbps;
B = zeros(Q, Q);
for i = 0 : Q - 1
    for j = 0 : Q - 1
        B(i + 1, j + 1) = sum(de2bi(bitxor(i, j)));
    end
end
