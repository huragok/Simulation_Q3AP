clear all;
close all;
clc;

addpath('./functions/')
%% 1. Generate the Gray mapped constellation
Nbps = 5;
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
    
    Eb2N0 = (0); % Eb/N0 in dB
    N = [200]; % Number of serial expansion
    n_Eb2N0 = length(Eb2N0);
    
    map_hans = cell(n_Eb2N0, 1);
    map_karim = cell(n_Eb2N0, 1);
    map_hans{1} = [30, 10, 12, 15, 8, 14, 31, 19, 23, 20, 22, 2, 25, 16, 24, 5, 13, 6, 29, 27, 21, 26, 28, 7, 0, 4, 1, 18, 17, 11, 9, 3;
                   30, 26, 12, 14, 13, 11, 27, 3, 21, 20, 22, 2, 29, 16, 25, 7, 24, 6, 23, 10, 9, 31, 28, 5, 0, 4, 1, 18, 17, 15, 8, 19] + 1;
    map_karim{1} = [11, 8, 9, 10, 23, 20, 21, 22, 3, 0, 1, 2, 15, 12, 13, 14, 27, 24, 25, 26, 7, 4, 5, 6, 19, 16, 17, 18, 31, 28, 29, 30;
                    11, 8, 9, 10, 23, 20, 21, 22, 3, 0, 1, 2, 15, 12, 13, 14, 27, 24, 25, 26, 7, 4, 5, 6, 19, 16, 17, 18, 31, 28, 29, 30] + 1; % A very simple remap based on "Trans-Modulation in Wireless Relay Networks"
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
BER_upperbound_hans = zeros(n_Eb2N0, 1);
BER_hans = zeros(n_Eb2N0, 1);
BER_upperbound_karim = zeros(n_Eb2N0, 1);
BER_karim = zeros(n_Eb2N0, 1);
matlabpool open 4
parfor i_Eb2N0 = 1 : n_Eb2N0
    tic;
    BER_upperbound_uniform(i_Eb2N0) = get_BER_upper_bound(X, map_uniform, mu_h, sigma2_h, sigma2_eps, sigma2_v(i_Eb2N0), N(i_Eb2N0), xi);
    BER_uniform(i_Eb2N0) = get_BER(X, map_uniform, mu_h, sigma2_h, sigma2_eps, sigma2_v(i_Eb2N0), M);
    BER_upperbound_hans(i_Eb2N0) = get_BER_upper_bound(X, map_hans{i_Eb2N0}, mu_h, sigma2_h, sigma2_eps, sigma2_v(i_Eb2N0), N(i_Eb2N0), xi);
    BER_hans(i_Eb2N0) = get_BER(X, map_hans{i_Eb2N0}, mu_h, sigma2_h, sigma2_eps, sigma2_v(i_Eb2N0), M);
    BER_upperbound_karim(i_Eb2N0) = get_BER_upper_bound(X, map_karim{i_Eb2N0}, mu_h, sigma2_h, sigma2_eps, sigma2_v(i_Eb2N0), N(i_Eb2N0), xi);
    BER_karim(i_Eb2N0) = get_BER(X, map_karim{i_Eb2N0}, mu_h, sigma2_h, sigma2_eps, sigma2_v(i_Eb2N0), M);
    
    disp(['Eb/N0 = ', num2str(Eb2N0(i_Eb2N0)), 'dB']);
    disp(['Uniform mapping, BER upperbound = ', num2str(BER_upperbound_uniform(i_Eb2N0))]);
    disp(['Uniform mapping, BER = ', num2str(BER_uniform(i_Eb2N0))]);
    disp(['Hans remapping, BER upperbound = ', num2str(BER_upperbound_hans(i_Eb2N0))]);
    disp(['Hans remapping, BER = ', num2str(BER_hans(i_Eb2N0))]);
    disp(['Karim remapping, BER upperbound = ', num2str(BER_upperbound_karim(i_Eb2N0))]);
    disp(['Karim remapping, BER = ', num2str(BER_karim(i_Eb2N0))]);
    toc;
    
    plot_mapping(X, map_hans{i_Eb2N0}(1, :), ['Eb/N0 = ', num2str(Eb2N0(i_Eb2N0)), 'dB, S-D']);
    plot_mapping(X, map_hans{i_Eb2N0}(2, :), ['Eb/N0 = ', num2str(Eb2N0(i_Eb2N0)), 'dB, R-D']);
end
matlabpool close

%% 4. Visualization
% range = 1 : 7;
% figure;
% semilogy(Eb2N0(range), BER_upperbound_uniform(range), 'bo--', 'linewidth', 2), hold on;
% semilogy(Eb2N0(range), BER_uniform(range), 'bo-', 'linewidth', 2), hold on;
% semilogy(Eb2N0(range), BER_upperbound_karim(range), 'rs--', 'linewidth', 2), hold on;
% semilogy(Eb2N0(range), BER_karim(range), 'rs-', 'linewidth', 2), hold on;
% semilogy(Eb2N0(range), BER_upperbound_hans(range), 'm^--', 'linewidth', 2), hold on;
% semilogy(Eb2N0(range), BER_hans(range), 'm^-', 'linewidth', 2), hold on;
% grid on;
% set(gca, 'Fontsize', 16);
% xlabel('E_b/N_0(dB)'), ylabel('BER'), legend('Uniform bound', 'Uniform sim.', 'Karim bound', 'Karim sim.', 'Hans bound', 'Hans sim.');