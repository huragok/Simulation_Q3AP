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
    
    map = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0;
           5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1; % A very simple remap based on "Trans-Modulation in Wireless Relay Networks"
    Eb2N0 = (-2 : 1 : 7); % Eb/N0 in dB
    N = [100, 200, 200, 300, 300, 400, 400, 800, 800, 1500]; % Number of serial expansion
    n_Eb2N0 = length(Eb2N0);
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
sigma2_v = pwr ./ (Nbps * 10 .^ (Eb2N0 / 10));

%% 3. Compute the BER
xi = 1 / 4;
M = 10 ^ 7;
map_uniform = [1 : Q; 1 : Q];

BER_upperbound_uniform = zeros(n_Eb2N0, 1);
BER_uniform = zeros(n_Eb2N0, 1);
BER_upperbound = zeros(n_Eb2N0, 1);
BER = zeros(n_Eb2N0, 1);
for i_Eb2N0 = 1 : n_Eb2N0
    tic;
    BER_upperbound_uniform(i_Eb2N0) = get_BER_upper_bound(X, map_uniform, mu_h, sigma2_h, sigma2_eps, sigma2_v(i_Eb2N0), N(i_Eb2N0), xi);
    BER_uniform(i_Eb2N0) = get_BER(X, map_uniform, mu_h, sigma2_h, sigma2_eps, sigma2_v(i_Eb2N0), M);
    BER_upperbound(i_Eb2N0) = get_BER_upper_bound(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v(i_Eb2N0), N(i_Eb2N0), xi);
    BER(i_Eb2N0) = get_BER(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v(i_Eb2N0), M);
    disp(['Eb/N0 = ', num2str(Eb2N0(i_Eb2N0)), 'dB']);
    disp(['Uniform mapping, BER upperbound = ', num2str(BER_upperbound_uniform(i_Eb2N0))]);
    disp(['Uniform mapping, BER = ', num2str(BER_uniform(i_Eb2N0))]);
    disp(['Simple remapping, BER upperbound = ', num2str(BER_upperbound(i_Eb2N0))]);
    disp(['Simple remapping, BER = ', num2str(BER(i_Eb2N0))]);
    toc;
end
%BER_upperbound = get_BER_upper_bound(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi)
%BER = get_BER(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v, M)

%% 4. Visualization
figure;
semilogy(Eb2N0, BER_upperbound_uniform, 'bo--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_uniform, 'bo-', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER_upperbound, 'rs--', 'linewidth', 2), hold on;
semilogy(Eb2N0, BER, 'rs-', 'linewidth', 2), hold on;
grid on;
set(gca, 'Fontsize', 16);
xlabel('E_b/N_0(dB)'), ylabel('BER'), legend('Uniform bound', 'Uniform sim.', 'Remap bound', 'Remap sim.');