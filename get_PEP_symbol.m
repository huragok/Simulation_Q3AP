function p = get_PEP_symbol(x, y, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi)
%   p = get_PEP_symbol(x, y, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi)
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

tau = tan(((1 : N)' - 1 / 2) * pi / N);
omega = xi + 1j * xi * tau;

mu0 = [mu_h(1); 0; 0];
R0 = diag([sigma2_h(1), sigma2_eps(1), sigma2_v]);
e0 = x(1) - y(1);
A0 = [abs(e0)^2, -e0' * y(1), e0';...
      -e0 * y(1)', abs(y(1))^2 - abs(x(1))^2, e0';...
      e0, e0, 0];

mu12 = [mu_h(2); mu_h(3); 0; 0; 0];
R12 = diag([sigma2_h(2), sigma2_h(3), sigma2_eps(2), sigma2_eps(3), sigma2_v]);
e1 = x(2) - y(2);
e2 = x(3) - y(3);
A12 = [abs(e1)^2, e1' * e2, -e1' * y(2), -e1' * y(3), e1';...
       e1 * e2', abs(e2)^2, -e2' * y(2), -e2' * y(3), e2';...
       -e1 * y(2)', -e2 * y(2)', abs(y(2))^2 - abs(x(2))^2, y(2)'*y(3) - x(2)'*x(3), e1';...
       -e1 * y(3)', -e2 * y(3)', y(2)*y(3)' - x(2)*x(3)', abs(y(3))^2 - abs(x(3))^2, e2';...
       e1, e2, e1, e2, 0];
   
psi = get_Psi_Gaussian(omega, mu0, R0, A0) .* get_Psi_Gaussian(omega, mu12, R12, A12);
p = sum(real(psi) + tau .* imag(psi)) / (2 * N);
end

function psi = get_Psi_Gaussian(omega, mu, R, A)
%   psi = get_Psi_Gaussian(omega, mu, R, A)
%   Compute the moment generating function (MGF) for the quadaratic form of
%   Gaussian random variables z'Az.
% _____________________________________________________________________________
%	Inputs:
% 		omega:      N-by-1 vector, a sequence of variables upon which the
%                   MGF is computed.
% 		
%       mu:         d-by-1 vector, the mean value of the Gaussian variable
%                   vector z
%       R:      d-by-d positive semi-definite matrix, the covariance
%                   matrix of the Gaussian variable vector z
%       A:          a d-by-d Hermitian matrix
%	Outputs:
%		psi:		N-by-1 vector, a sequence of MGF values evaluated
% _____________________________________________________________________________

d = length(mu); % Size of the Gaussian random vector
psi = zeros(N, 1);
for n = 1 : N
    IRA = eye(d) + omega(n) * R * A;
    psi(n) = exp(-omega(n) * mu' * A / IRA * mu) / det(IRA);
end

end
