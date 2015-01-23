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
    
    Eb2N0 = (-2 : 1 : 7); % Eb/N0 in dB
    N = [100, 200, 200, 300, 300, 400, 400, 800, 800, 1500]; % Number of serial expansion
    n_Eb2N0 = length(Eb2N0);
    
    map_hans = cell(n_Eb2N0, 1);
    map_karim = cell(n_Eb2N0, 1);
    map_hans{1} = [11, 10, 7, 9, 15, 2, 4, 8, 3, 13, 0, 1, 6, 14, 12, 5;
                   11, 5, 8, 9, 14, 6, 4, 13, 3, 2, 15, 1, 7, 10, 12, 0] + 1;
    map_hans{2} = [10, 9, 7, 11, 14, 15, 12, 8, 6, 3, 4, 5, 13, 1, 0, 2;
                   5, 9, 8, 11, 14, 10, 12, 13, 6, 7, 4, 1, 2, 0, 15, 3] + 1;
    map_hans{3} = [14, 11, 12, 13, 5, 9, 8, 10, 2, 1, 15, 3, 6, 7, 4, 0;
                   14, 15, 12, 9, 10, 8, 7, 11, 13, 1, 0, 3, 6, 2, 4, 5] + 1;
    map_hans{4} = [10, 11, 9, 7, 6, 15, 12, 4, 13, 3, 1, 0, 14, 2, 5, 8;
                   5, 11, 9, 8, 6, 10, 13, 0, 2, 3, 1, 15, 14, 7, 4, 12] + 1;
    map_karim{1} = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0;
                    5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1; % A very simple remap based on "Trans-Modulation in Wireless Relay Networks"
    map_karim{2} = map_karim{1};
    map_karim{3} = map_karim{1};
    map_karim{4} = map_karim{1};
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
BER_upperbound_hans = zeros(n_Eb2N0, 1);
BER_hans = zeros(n_Eb2N0, 1);
BER_upperbound_karim = zeros(n_Eb2N0, 1);
BER_karim = zeros(n_Eb2N0, 1);
for i_Eb2N0 = 1 : 4
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
    plot_mapping(X, map_hans{i_Eb2N0}(1, :), ['Eb/N0 = ', num2str(Eb2N0(i_Eb2N0)), 'dB, R-D']);
end

%% 4. Visualization
range = 1 : 4;
figure;
semilogy(Eb2N0(range), BER_upperbound_uniform(range), 'bo--', 'linewidth', 2), hold on;
semilogy(Eb2N0(range), BER_uniform(range), 'bo-', 'linewidth', 2), hold on;
semilogy(Eb2N0(range), BER_upperbound_karim(range), 'rs--', 'linewidth', 2), hold on;
semilogy(Eb2N0(range), BER_karim(range), 'rs-', 'linewidth', 2), hold on;
semilogy(Eb2N0(range), BER_upperbound_hans(range), 'm^--', 'linewidth', 2), hold on;
semilogy(Eb2N0(range), BER_hans(range), 'm^-', 'linewidth', 2), hold on;
grid on;
set(gca, 'Fontsize', 16);
xlabel('E_b/N_0(dB)'), ylabel('BER'), legend('Uniform bound', 'Uniform sim.', 'Karim bound', 'Karim sim.', 'Hans bound', 'Hans sim.');