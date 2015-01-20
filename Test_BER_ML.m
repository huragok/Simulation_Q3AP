clear all;
close all;
clc;

%% 1. Generate the Gray mapped constellation
Nbps = 4;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);
Q = 2 ^ Nbps;

%% 2. Channel settings (Let us try AWGN first)
channel = 'AWGN'; % Channel model, can be specified as AWGN, Rayleigh, Rician, Rayleigh_imp, Rician_imp
if strcmp(channel, 'AWGN') % AWGN channel
    mu_h = [1; 1; 1];
	sigma2_h = [0; 0; 0];
	sigma2_eps = [0; 0; 0];
    
    %map = [0, 12, 1, 13, 6, 5, 2, 14, 3, 10, 8, 9, 4, 7, 15, 11;
           %0, 6, 1, 5, 9, 13, 7, 11, 15, 10, 12, 3, 14, 2, 8, 4] + 1; % This is the initial result for the AWGN channel at 5dB, probably wrong due to the ordering issue
elseif strcmp(channel, 'Rayleigh') % Rayleigh fading channel with perfect CSIR
    mu_h = [0; 0; 0];
    sigma2_h = [1; 1; 1];
    sigma2_eps = [0; 0; 0];
elseif strcmp(channel, 'Rician') % Rician fading channel with perfect CSIR
    K = 10;
    mu_h = sqrt(K / (K + 1)) * ones(3, 1);
    sigma2_h = 1 / (K + 1) * ones(3, 1);
    sigma2_eps = [0; 0; 0];
elseif strcmp(channel, 'Rayleigh_imp')% Rayleigh fading channel with imperfect CSIR
    mu_h = [0; 0; 0];
    sigma2_h = [1; 1; 1];
    sigma2_eps = [0.2; 0.2; 0.2];
elseif strcmp(channel, 'Rician_imp') % Rician fading channel with perfect CSIR
    K = 10;
    mu_h = sqrt(K / (K + 1)) * ones(3, 1);
    sigma2_h = 1 / (K + 1) * ones(3, 1);
    sigma2_eps = [0.2; 0.2; 0.2];
else
    error('Wrong channel specified!')
end
Eb2N0 = 0; % Eb/N0 in dB
sigma2_v = pwr / (Nbps * 10 ^ (Eb2N0 / 10));

%% 3. Compute the BER
map = [1 : Q; 1 : Q];
M = 10000; % scalar, the number of ML decoding carried out

%% 4. Generate the channels and the received symbols for the demodulator
constellation = X([1 : Q; map]); % The constellation vectors (for all 3 transmissions)
i_symbols_mod = unidrnd(Q, [1, M]) - 1; % The symbol index (0 ~ Q-1) of the actually transmitted symbols
s = constellation(:, i_symbols_mod + 1); % The transmitted symbols

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

%% 5. Run the ML demodulator
d = zeros(Q, M);
for q = 1 : Q
    r = zeros(2, M);
    r(1, :) = h_est(1, :) * constellation(1, q);
    r(2, :) = sum(h_est(2 : 3, :) .* (constellation(2 : 3, q) * ones(1, M)));
    d(q, :) = sum((y - r) .^ 2);
end
[~, i_symbols_demod] = min(d);
i_symbols_demod = i_symbols_demod - 1;

%% Count BER
B = get_n_diff_bits(Nbps);
n_diff_bits = zeros(1, M);
for m = 1 : M
    n_diff_bits(m) = B(i_symbols_mod(m) + 1, i_symbols_demod(m) + 1);
end
BER = sum(n_diff_bits) / (M * Nbps)