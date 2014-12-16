function delta = rand_delta(x, y, mu_h, sigma2_h, sigma2_eps, sigma2_v, M)
%   delta = rand_delta(x, y, mu_h, sigma2_h, sigma2_eps, sigma2_v, M)
%   Random variable generator for delta. Used to verify the approximated
%   PEP of 2 symbols.
% _____________________________________________________________________________
%	Inputs:
% 		x:          3-by-1 vector, the actually transmitted symbols by the 3
%                   different mappings.
%       y:          3-by-1 vector, the alternating transmitted symbols by
%                   the 3 different mappings.
%       mu_h:       3-by-1 vector, the mean value of the Rician channels
%                   (LOS component)
%       sigma2_h:   3-by-1 vector, the variance of the Rician channels
%                   (fading component)
%       sigma2_eps: 3-by-1 vector, the variance of the channel estimation
%                   error
%       sigma2_v:   scalar, the variance of the received AWGN noise at the
%                   destination
%       M:          scalar, the number of random delta to generate
%	Outputs:
%		delta:      M-by-1 vector, the randomly generated delta
% _____________________________________________________________________________
% Author: Wenhao Wu
% Email: wnhwu@ucdavis.edu
% Date: 12/16/2014
% Codename: Dunkirk
% _____________________________________________________________________________
% References:
%   [1] Harvind Singh Samra. Effective signal processing algorithms for packet
%       re-transmission diversity. Ph.D. thesis, Dept. Elect. Eng., UC Davis, 
%       CA, 2004.
% _____________________________________________________________________________
h0 = 

end

function r = crandn(mu, sigma2, M)
%   r = crandn(mu, sigma2, M)
%   Complex circularly symmetric Gaussian random variable generator.
% _____________________________________________________________________________
%	Inputs:
% 		mu:     scalar, the mean value of the complex Gaussian variable
%       sigma2: scalar, the covariance of the complex Gaussian variable,
%               equally divided between the real and the image part
%       M:      scalar, the number of random delta to generate
%   Outputs:
%		r:      M-by-1 vector, the randomly generated complex Gaussian
%               random variables
% _____________________________________________________________________________


end