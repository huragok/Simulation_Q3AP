function BER = get_BER_upper_bound(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi)
%   BER = get_BER_upper_bound(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi)
%   Compute the BER upper bound given the mapping scheme and the channel
%   conditions
% _________________________________________________________________________
%	Inputs:
%       X:          Q-by-1 vector, the constellation symbols.
% 		map:        2-by-Q matrix, each row is a permutation of 1 : Q
%                   indicating mapping of the symbols in the retransmission
%                   round from the source and the relay
%       mu_h:       3-by-1 vector, the mean value of the Rician channels
%                   (LOS component)
%       sigma2_h:   3-by-1 vector, the variance of the Rician channels
%                   (fading component)
%       sigma2_eps: 3-by-1 vector, the variance of the channel estimation
%                   error
%       sigma2_v:   scalar, the variance of the received AWGN noise at the
%                   destination
%       N:          scalar, the number of integration points when
%                   approximating the PEP that controls the accuracy
%       xi:         scalar, parameters in the numerical integration that
%                   must ensures convergence. 1/4 is recommended
%	Outputs:
%		BER:		scalar, the upper bound of BER computed with pairwise
%                   error probability 
% _________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 01/14/2015
% Codename: Dunkirk
% _________________________________________________________________________
% References:
%   [1] Harvind Singh Samra. Effective signal processing algorithms for 
%       packet re-transmission diversity. Ph.D. thesis, Dept. Elect. Eng.,
%       UC Davis, CA, 2004.
% _________________________________________________________________________

Q = length(X);
Nbps = log2(Q); 
B = get_n_diff_bits(Nbps);

g = zeros(Q); %g[p, q] is a pair wise error probability (in terms of bits) between symbol p and q
for p = 1 : Q
    x = X([p; map(:, p)]); % The actually transmitted 3 symbols
    for q = 1 : Q
        if p ~= q
            y = X([q; map(:, q)]); % The alternative 3 symbols
            g(p, q) = B(p, q) * get_PEP_symbol(x, y, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi) / Q;
        end
    end
end

BER = sum(sum(g));
