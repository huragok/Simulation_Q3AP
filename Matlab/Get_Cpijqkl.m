clear all;
close all;
clc;

addpath('./functions/');
%% 1. Generate the Gray mapped constellation
Nbps = 2;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);

%% 2. Generate all the Q ^ 6 index vectors [p, i, j, q, k, l]' (or equavilently [p, a, c, q, b, d]')
Q = 2 ^ Nbps;
order = [4, 1, 5, 2, 6, 3];

%% 3. Channel settings (Let us try AWGN first)
channel = 'Rician'; % Channel model, can be specified as AWGN, Rayleigh, Rician, Rayleigh_imp, Rician_imp
if strcmp(channel, 'AWGN') % AWGN channel
    mu_h = [1; 1; 1];
	sigma2_h = [0; 0; 0];
	sigma2_eps = [0; 0; 0];
    
    Eb2N0 = 0; % Eb/N0 in dB
    N = 200; % Number of serial expansion
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
    Eb2N0 = 0; % Eb/N0 in dB
    N = 200; % Number of serial expansion
    n_Eb2N0 = length(Eb2N0);
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

%% 4. Compute PEP for the symbols
xi = 1 / 4;

%idxs_cell = num2cell(idxs, 1);
B = get_n_diff_bits(Nbps);

tic;
matlabpool open 4 % My computers has a Core i7-4790 CPU with 4 cores and it takes ~8000 sec. If your computer has more cores it is possible to open more thread to further speed it up
for i_Eb2N0 = 1 : 1
    tic;
    sigma2_v_tmp = sigma2_v(i_Eb2N0);
    N_tmp = N(i_Eb2N0);
    
    PEP_MGF_bit = zeros(1, Q ^ 6); % PEP computed with MGF method
    for q = 1 : Q ^ 6
        idxs = get_idxs(q, Q, order);
        symbols_base = X(idxs(1 : 3)); % The actually transmitted 3 symbols
        symbols_alt = X(idxs(4 : 6)); % The alternative 3 symbols
        if all(idxs(1 : 3) ~= idxs(4 : 6))
            PEP_MGF = get_PEP_symbol(symbols_base, symbols_alt, mu_h, sigma2_h, sigma2_eps, sigma2_v_tmp, N_tmp, xi);
            PEP_MGF_bit(q) =  PEP_MGF * B(idxs(1), idxs(4)) / Q;
        end
    end
    
    filename = ['q3aptest_', num2str(Q), type_mod, '_', num2str(Eb2N0(i_Eb2N0)), 'dB.data'];
    fileID = fopen(filename, 'w+');
    fprintf(fileID, '  %18.16e', PEP_MGF_bit);
    toc;
    fclose(fileID);
    
    disp(['Eb/N0 = ', num2str(Eb2N0(i_Eb2N0)), 'dB completed']);
end
matlabpool close

