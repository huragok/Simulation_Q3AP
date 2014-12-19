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
for q = 0 : Q ^ 6 - 1
    q_residual = q;
    for d = 1 : 6
        idxs(d, q + 1) = mod(q_residual, Q);
        q_residual = floor(q_residual / Q);
    end
end
idxs = idxs + 1;

symbols_base = num2cell(X(idxs([1, 3, 5], :)), 1); % The actually transmitted 3 symbols
symbols_alt = num2cell(X(idxs([2, 4, 6], :)), 1); % The alternative 3 symbols

% Randomly pick a pair of symbols
% symbol_base = symbols_base(:, 1);
% symbol_alt = symbols_alt(:, 1);

%% 3. Channel settings (Let us try AWGN first)
% AWGN channel
% mu_h = [1; 1; 1];
% sigma2_h = [0; 0; 0];
% sigma2_eps = [0; 0; 0];

% Rayleigh fading channel with perfect CSIR
% mu_h = [0; 0; 0];
% sigma2_h = [1; 1; 1];
% sigma2_eps = [0; 0; 0];

% Rician fading channel with perfect CSIR
% K = 10;
% mu_h = sqrt(K / (K + 1)) * ones(3, 1);
% sigma2_h = 1 / (K + 1) * ones(3, 1);
% sigma2_eps = [0; 0; 0];

% Rayleigh fading channel with imperfect CSIR
% mu_h = [0; 0; 0];
% sigma2_h = [1; 1; 1];
% sigma2_eps = [0.2; 0.2; 0.2];

% Rician fading channel with perfect CSIR
K = 10;
mu_h = sqrt(K / (K + 1)) * ones(3, 1);
sigma2_h = 1 / (K + 1) * ones(3, 1);
sigma2_eps = [0.2; 0.2; 0.2];

Eb2N0 = 5; % Eb/N0 in dB
sigma2_v = pwr / Nbps * 10 ^ (Eb2N0 / 10);

%% 4. Compute PEP
N = 64;
xi = 1 / 4;
M = 50000;
tic;

PEP_MGF = NaN(1, Q ^ 6); % PEP computed with MGF method
PEP_MC = NaN(1, Q ^ 6); % PEP computed with MGF method

idxs_cell = num2cell(idxs, 1);
matlabpool open 4
parfor q = 1 : Q ^ 6
    if all(idxs_cell{q}([1, 3, 5]) ~= idxs_cell{q}([2, 4, 6]))
        PEP_MGF(q) = get_PEP_symbol(symbols_base{q}, symbols_alt{q}, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi);
        PEP_MC(q) = get_PEP_symbol_MC(symbols_base{q}, symbols_alt{q}, mu_h, sigma2_h, sigma2_eps, sigma2_v, M);
    end
end
matlabpool close

toc;

%% 5. Visualization
[PEP_MGF, I] = sort(PEP_MGF);
PEP_MC = PEP_MC(I);

PEP_MFG_sorted = PEP_MGF(~isnan(PEP_MGF));
PEP_MC_sorted = PEP_MC(~isnan(PEP_MC));

l = length(PEP_MFG_sorted);

figure;
plot(1 : l, PEP_MFG_sorted, 'b-', 'lineWidth', 2), hold on;
plot(1 : l, PEP_MFG_sorted - PEP_MC_sorted, 'r:', 'lineWidth', 2), hold on;
set(gca, 'fontsize', 16), xlabel('Index'), ylabel('PEP'), grid on;
legend('MGF', 'MGF - FC')


