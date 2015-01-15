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
N = [50 : 10 : 180]; % Compute the BER upperbound at different numbers of serial expansion
xi = [1 / 4];
BER_upperbound = cell(length(xi), 1);
%matlabpool open 3
for i_xi = 1 : length(xi)
    BER_upperbound{i_xi} = zeros(length(N), 1);
    for i_N = 1 : length(N)
       BER_upperbound{i_xi}(i_N) = get_BER_upper_bound(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v, N(i_N), xi(i_xi));
       disp(['xi = ', num2str(xi(i_xi)), ', N = ', num2str(N(i_N)), ' completed.']);
    end
end
%matlabpool close

%% Visualization
figure;
cmap = colormap(hsv(length(xi)));
legend_item = cell(length(xi), 1);
for i_xi = 1 : length(xi)
    plot(N, BER_upperbound{i_xi}, 'linewidth', 2, 'Color', cmap(i_xi, :)), hold on;
    legend_item{i_xi} = ['\xi = ', num2str(xi(i_xi))];
end
legend(legend_item);
grid on;
set(gca, 'fontsize', 18), xlabel('N'), ylabel('BER upperbound')