function get_BER_upper_bound_from_C(map, C)
%   BER = get_BER_upper_bound_from_C(map, C)
%   Compute the BER upper bound given the mapping scheme and the cost
%   matrix
% _________________________________________________________________________
%	Inputs:
% 		map:        2-by-Q matrix, each row is a permutation of 1 : Q
%                   indicating mapping of the symbols in the retransmission
%                   round from the source and the relay
%       C:          Q ^ 6-by-1 cost matrix. In the order of C(p, u(p), v(p)
%                   , q, u(q), v(q)). Right most order increases first.
%	Outputs:
%		BER:		scalar, the upper bound of BER computed with pairwise
%                   error probability 
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 01/23/2015
% Codename: Dunkirk
% _________________________________________________________________________
% References:
%   [1] Harvind Singh Samra. Effective signal processing algorithms for 
%       packet re-transmission diversity. Ph.D. thesis, Dept. Elect. Eng.,
%       UC Davis, CA, 2004.
% _________________________________________________________________________

Q = size(map, 2);
map_complete = [1 : Q; map] - 1; % The complete mapping of all the 3 transmissions, each row is a permutation of 0 : Q - 1

BER = 0;
for p = 1 : q
    for q = 1 : q
        % Compute the index corresponding to the 6-D coordinate (p, u(p), v(p), q, u(q), v(q))
        idx = Q ^ 6 - Q .^ (5 : -1 : 0) * [map_complete(:, p); map_complete(:, q)];
        BER = BER + C(idx);
    end
end
