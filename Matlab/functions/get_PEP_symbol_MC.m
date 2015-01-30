function p = get_PEP_symbol_MC(x, y, mu_h, sigma2_h, sigma2_eps, sigma2_v, M)
%   p = get_PEP_symbol_MC(x, y, mu_h, sigma2_h, sigma2_eps, sigma2_v, M)
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
%       M:          scalar, the number of random delta used in Monte-Carlo
%                   simulation
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

delta = rand_delta(x, y, mu_h, sigma2_h, sigma2_eps, sigma2_v, M);
p = sum(delta < 0) / M;

end

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

h0 = crandn(mu_h(1), sigma2_h(1), M);
eps0 = crandn(0, sigma2_eps(1), M);
v0 = crandn(0, sigma2_v, M);

h1 = crandn(mu_h(2), sigma2_h(2), M);
h2 = crandn(mu_h(3), sigma2_h(3), M);
eps1 = crandn(0, sigma2_eps(2), M);
eps2 = crandn(0, sigma2_eps(3), M);
v1 = crandn(0, sigma2_v, M);

delta = abs((x(1) - y(1)) * h0 - y(1) * eps0 + v0) .^ 2 ...
      - abs(-x(1) * eps0 + v0) .^ 2 ...
      + abs((x(2) - y(2)) * h1 + (x(3) - y(3)) * h2 - y(2) * eps1 - y(3) * eps2 + v1) .^ 2 ...
      - abs(-x(2) * eps1 - x(3) * eps2 + v1) .^ 2;
end
