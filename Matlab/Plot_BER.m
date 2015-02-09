clear all;
close all;
clc;

addpath('./functions/')
%% 1. Generate the Gray mapped constellation
Nbps = 4;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);
Q = 2 ^ Nbps;

%% 2. Channel settings (Let us try AWGN first)
%channel = 'AWGN'; % Channel model, can be specified as AWGN, Rayleigh, Rician, Rayleigh_imp, Rician_imp
% channel = 'Rayleigh';
channel = 'Rician';
if strcmp(channel, 'AWGN') % AWGN channel
    mu_h = [1; 1; 1];
	sigma2_h = [0; 0; 0];
	sigma2_eps = [0; 0; 0];
    
    Eb2N0 = (-2 : 1 : 7); % Eb/N0 in dB
    N = [100, 200, 200, 300, 400, 400, 400, 800, 800, 1500]; % Number of serial expansion
    n_Eb2N0 = length(Eb2N0);
    
    map_hans = cell(n_Eb2N0, 1);
    map_karim = cell(n_Eb2N0, 1);
    map_hans{1} = [13, 11, 7, 3, 4, 8, 12, 0, 15, 9, 5, 1, 14, 10, 6, 2;
                   7, 9, 13, 1, 12, 8, 4, 0, 5, 11, 15, 3, 6, 10, 14, 2] + 1;
    map_hans{2} = [5, 1, 15, 11, 6, 2, 14, 10, 7, 3, 13, 9, 4, 0, 12, 8;
                   15, 3, 5, 9, 14, 2, 6, 10, 13, 1, 7, 11, 12, 0, 4, 8] + 1;
    map_hans{3} = [15, 12, 5, 6, 11, 8, 9, 10, 7, 4, 13, 14, 3, 0, 1, 2;
                   5, 4, 15, 14, 9, 8, 11, 10, 13, 12, 7, 6, 1, 0, 3, 2] + 1;
    map_hans{4} = [13, 12, 5, 6, 11, 8, 9, 10, 7, 4, 15, 14, 3, 0, 1, 2;
                   7, 4, 15, 14, 9, 8, 11, 10, 13, 12, 5, 6, 1, 0, 3, 2] + 1;
    map_hans{5} = [15, 11, 13, 3, 6, 2, 12, 8, 5, 9, 7, 1, 4, 0, 14, 10;
                   5, 9, 7, 1, 14, 2, 4, 8, 15, 11, 13, 3, 12, 0, 6, 10] + 1;
    map_karim{1} = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0;
                    5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1; % A very simple remap based on "Trans-Modulation in Wireless Relay Networks"
    map_karim{2} = map_karim{1};
    map_karim{3} = map_karim{1};
    map_karim{4} = map_karim{1};
    map_karim{5} = map_karim{1};
elseif strcmp(channel, 'Rayleigh') % Rayleigh fading channel with perfect CSIR
    mu_h = [0; 0; 0];
    sigma2_h = [1; 1; 1];
    sigma2_eps = [0; 0; 0];
    
    Eb2N0 = (-2 : 1 : 4); % Eb/N0 in dB
    N = [100, 200, 200, 300, 400, 400, 400]; % Number of serial expansion
    n_Eb2N0 = length(Eb2N0);
    map_hans = cell(n_Eb2N0, 1);
    map_karim = cell(n_Eb2N0, 1);
    
    map_hans{1} = [15, 8, 14, 11, 2, 5, 3, 6, 12, 9, 13, 10, 1, 4, 0, 7;
                   0, 7, 12, 11, 5, 2, 1, 6, 4, 3, 8, 15, 9, 14, 13, 10] + 1;
    map_hans{2} = [15, 2, 13, 0, 8, 1, 9, 12, 7, 3, 5, 4, 10, 6, 14, 11;
                   4, 5, 0, 7, 13, 10, 9, 14, 8, 1, 12, 3, 15, 6, 11, 2] + 1;
    map_hans{3} = [0, 2, 5, 14, 12, 13, 4, 8, 1, 3, 7, 15, 11, 6, 9, 10;
                   15, 7, 9, 11, 5, 2, 1, 14, 13, 6, 8, 10, 12, 0, 4, 3] + 1;
    map_hans{4} = [8, 13, 11, 6, 5, 0, 12, 1, 9, 4, 10, 15, 14, 3, 7, 2;
                   7, 0, 15, 3, 10, 12, 11, 6, 13, 1, 5, 2, 9, 4, 8, 14] + 1;
    map_hans{5} = [7, 0, 15, 3, 10, 12, 11, 6, 13, 1, 5, 2, 9, 4, 8, 14;
                   0, 5, 3, 14, 13, 8, 4, 9, 1, 12, 2, 7, 6, 11, 15, 10] + 1;
    map_hans{6} = [14, 3, 10, 13, 9, 4, 15, 8, 2, 5, 6, 11, 7, 0, 1, 12;
                   5, 4, 15, 8, 6, 3, 2, 1, 7, 0, 13, 12, 10, 9, 14, 11] + 1;
    map_hans{7} = [11, 4, 15, 6, 1, 0, 5, 2, 14, 13, 7, 3, 9, 12, 8, 10;
                   5, 4, 8, 1, 9, 7, 10, 14, 13, 0, 12, 3, 11, 2, 15, 6] + 1;
               
    map_karim{1} = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0;
                    5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1; % A very simple remap based on "Trans-Modulation in Wireless Relay Networks"
    map_karim{2} = map_karim{1};
    map_karim{3} = map_karim{1};
    map_karim{4} = map_karim{1};
    map_karim{5} = map_karim{1};
    map_karim{6} = map_karim{1};
    map_karim{7} = map_karim{1};
elseif strcmp(channel, 'Rician') % Rician fading channel with perfect CSIR
    K = 10;
    mu_h = sqrt(K / (K + 1)) * ones(3, 1);
    sigma2_h = 1 / (K + 1) * ones(3, 1);
    sigma2_eps = [0; 0; 0];
    
    Eb2N0 = (-2 : 1 : 6); % Eb/N0 in dB
    N = [100, 200, 200, 300, 400, 400, 400, 800, 800]; % Number of serial expansion
    n_Eb2N0 = length(Eb2N0);
    map_hans = cell(n_Eb2N0, 1);
    map_karim = cell(n_Eb2N0, 1);
    
    map_hans{1} = [15, 11, 7, 3, 12, 8, 4, 0, 13, 9, 5, 1, 14, 10, 6, 2;
                   5, 9, 13, 1, 4, 8, 12, 0, 7, 11, 15, 3, 6, 10, 14, 2] + 1;
    map_hans{2} = [7, 3, 15, 11, 6, 2, 14, 10, 5, 1, 13, 9, 4, 0, 12, 8;
                   13, 1, 5, 9, 14, 2, 6, 10, 15, 3, 7, 11, 12, 0, 4, 8] + 1;
    map_hans{3} = [5, 4, 13, 6, 1, 0, 3, 2, 7, 12, 15, 14, 9, 8, 11, 10;
                   15, 12, 13, 14, 3, 0, 1, 2, 7, 4, 5, 6, 11, 8, 9, 10] + 1;
    map_hans{4} = [7, 4, 5, 6, 11, 8, 9, 10, 15, 12, 13, 14, 3, 0, 1, 2
                   13, 12, 5, 14, 9, 8, 11, 10, 15, 4, 7, 6, 1, 0, 3, 2] + 1;
    map_hans{5} = [5, 1, 7, 9, 4, 0, 12, 8, 13, 3, 15, 11, 6, 2, 14, 10;
                   15, 3, 7, 11, 12, 0, 4, 8, 13, 1, 5, 9, 14, 2, 6, 10] + 1;
    map_hans{6} = [5, 14, 7, 12, 9, 2, 11, 0, 15, 4, 13, 6, 1, 8, 3, 10
                   5, 14, 7, 12, 11, 2, 9, 0, 15, 4, 13, 6, 3, 8, 1, 10] + 1;
    map_hans{7} = [7, 11, 13, 1, 12, 2, 6, 8, 5, 9, 15, 3, 14, 0, 4, 10;
                   7, 11, 13, 1, 12, 2, 6, 8, 5, 9, 15, 3, 14, 0, 4, 10] + 1;
    map_hans{8} = [5, 9, 15, 3, 14, 0, 4, 10, 7, 11, 13, 1, 12, 2, 6, 8;
                   5, 9, 15, 3, 14, 0, 4, 10, 7, 11, 13, 1, 12, 2, 6, 8] + 1;
    map_hans{9} = [13, 3, 5, 11, 6, 8, 14, 0, 15, 1, 7, 9, 4, 10, 12, 2;
                   13, 3, 5, 11, 6, 8, 14, 0, 15, 1, 7, 9, 4, 10, 12, 2] + 1;
    
    map_karim{1} = [5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0;
                    5, 6, 7, 4, 9, 10, 11, 8, 13, 14, 15, 12, 1, 2, 3, 0] + 1; % A very simple remap based on "Trans-Modulation in Wireless Relay Networks"
    map_karim{2} = map_karim{1};
    map_karim{3} = map_karim{1};
    map_karim{4} = map_karim{1};
    map_karim{5} = map_karim{1};
    map_karim{6} = map_karim{1};
    map_karim{7} = map_karim{1};
    map_karim{8} = map_karim{1};
    map_karim{9} = map_karim{1};
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
range = 1 : 7;
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