clear all;
close all;
clc;

%% 1. Generate the Gray mapped constellation
Nbps = 2;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);

%% 2. Generate all the Q ^ 6 index vectors [p, q, i, k, j, n]' (or equavilently [p, q, a, b, c, d]')
Q = 2 ^ Nbps;
idxs = zeros(6, Q ^ 6);
for q = 0 : Q ^ 6 - 1;
    q_residual = q;
    for d = 1 : 6
        idxs(d, q + 1) = mod(q_residual, Q);
        q_residual = floor(q_residual / Q);
    end
end
idxs = idxs + 1;

symbols_base = X(idxs([1, 3, 5], :)); % The actually transmitted 3 symbols
symbols_alt = X(idxs([2, 4, 6], :)); % The alternative 3 symbols

% Randomly pick a pair of symbols
symbol_base = symbols_base(:, 1);
symbol_alt = symbols_alt(:, 1);

%% 3. Channel settings (Let us try AWGN first)
mu_h = [1; 1; 1];
sigma2_h = [0; 0; 0];
sigma2_eps = [0; 0; 0];
Eb2N0 = 5; % Eb/N0 in dB
sigma2_v = pwr / Nbps * 10 ^ (Eb2N0 / 10);

%% 4. Try MGF approximation method
N = 64;
xi = 1 / 4;
tic;
p = get_PEP_symbol(symbol_base, symbol_alt, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi);
toc;
disp(['P{delta < 0} = ', num2str(p)])

%% 5. Try Monte-Carlo Simulation method
M = 10000;
tic;
p_MC = get_PEP_symbol_MC(symbol_base, symbol_alt, mu_h, sigma2_h, sigma2_eps, sigma2_v, M);
toc;
disp(['P{delta < 0} = ', num2str(p_MC), ' (Monte-Carlo method)'])
