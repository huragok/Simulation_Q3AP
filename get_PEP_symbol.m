function p = get_PEP_symbol(x, y, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi)
%   p = get_PEP_symbol(x, y, mu_h, sigma2_h, sigma2_eps, sigma2_v)
%   Compute the pairwise error probability between two set of symbols in
%   the realy - HARQ model with Rician channel and imperfect CSIR
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
%       N:          scalar, the number of integration points when
%                   approximating the PEP that controls the accuracy
%       xi:         scalar, parameters in the numerical integration that
%                   must ensures convergence. 1/4 is recommended
%	Outputs:
%		p:			scalar, the pairwise error probability of symbols in
%                   relay - HARQ network

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

end

function psi = get_Psi_Gaussian(omega, mu, Sigma, A)
%   psi = get_Psi_Gaussian(omega, mu, Sigma, A)
%   Compute the moment generating function (MGF) for the quadaratic form of
%   Gaussian random variables z'Az.
% _____________________________________________________________________________
%	Inputs:
% 		omega:      N-by-1 vector, a sequence of variables upon which the
%                   MGF is computed.
% 		
%       mu:         d-by-1 vector, the mean value of the Gaussian variable
%                   vector z
%       Sigma:      d-by-d positive semi-definite matrix, the covariance
%                   matrix of the Gaussian variable vector z
%       A:          a d-by-d Hermitian matrix
%	Outputs:
%		psi:		N-by-1 vector, a sequence of MGF values evaluated
% _____________________________________________________________________________

end
