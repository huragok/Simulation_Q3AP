function p = get_BER_upper_bound(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi)
%   p = get_BER_upper_bound(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi)
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
%		p:			scalar, the upper bound of BER computed with pairwise
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