function p = get_BER(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v, M)
%   p = get_BER(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v, M)
%   Compute the BER by actually running a ML decoder M times given the 
%   mapping scheme and the channel conditions.
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
%       M:          scalar, the number of ML decoding carried out
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