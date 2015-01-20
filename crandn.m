function r = crandn(mu, sigma2, M, N)
%   r = crandn(mu, sigma2, M)
%   Complex circularly symmetric Gaussian random variable generator.
% _____________________________________________________________________________
%	Inputs:
% 		mu:     scalar, the mean value of the complex Gaussian variable
%       sigma2: scalar, the covariance of the complex Gaussian variable,
%               equally divided between the real and the image part
%       M, N:   scalar, the size of random delta to generate
%   Outputs:
%		r:      M-by-1 vector, the randomly generated complex Gaussian
%               random variables
% _____________________________________________________________________________

r = mu + sqrt(sigma2 / 2) * (randn(M, N) + 1j * randn(M, N));
end