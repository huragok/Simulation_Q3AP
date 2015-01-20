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
%		p:			scalar, the BER computed with Monte-Carlo simulation by
%                   running a ML demodulator 
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
Nbps = round(log2(Q));
constellation = X([1 : Q; map]); % The constellation vectors (for all 3 transmissions)
i_symbols_mod = unidrnd(Q, [1, M]); % The symbol index (1 ~ Q) of the actually transmitted symbols
s = constellation(:, i_symbols_mod); % The transmitted symbols

h = zeros(3, M); % The actual channels
for i_transmission = 1 : 3
    h(i_transmission, :) = crandn(mu_h(i_transmission), sigma2_h(i_transmission), 1, M);
end

delta_h = zeros(3, M); % The channel estimation errors
for i_transmission = 1 : 3
    delta_h(i_transmission, :) = crandn(0, sigma2_eps(i_transmission), 1, M);
end
h_est = h + delta_h; % The estimated channels

y = zeros(2, M);
y(1, :) = h(1, :) .* s(1, :);
y(2, :) = sum(h(2 : 3, :) .* s(2 : 3, :), 1);
y = y + crandn(0, sigma2_v, 2, M); % The received signals

d = zeros(Q, M);
for q = 1 : Q
    r = zeros(2, M);
    r(1, :) = h_est(1, :) * constellation(1, q);
    r(2, :) = sum(h_est(2 : 3, :) .* (constellation(2 : 3, q) * ones(1, M)));
    d(q, :) = sum(abs(y - r) .^ 2);
end
[~, i_symbols_demod] = min(d);

B = get_n_diff_bits(Nbps);
n_diff_bits = zeros(1, M);
for m = 1 : M
    n_diff_bits(m) = B(i_symbols_mod(m), i_symbols_demod(m));
end
p = sum(n_diff_bits) / (M * Nbps);