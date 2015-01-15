clear all;
close all;
clc;

%% 1. Generate the Gray mapped constellation
Nbps = 4;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);

%% 2. Generate all the Q ^ 6 index vectors [p, i, j, q, k, l]' (or equavilently [p, a, c, q, b, d]')
Q = 2 ^ Nbps;
idxs = zeros(6, Q ^ 6);
for q = 0 : Q ^ 6 - 1
    q_residual = q;
    for d = 6 : -1 : 1
        idxs(d, q + 1) = mod(q_residual, Q);
        q_residual = floor(q_residual / Q);
    end
end
idxs = idxs + 1;

symbols_base = num2cell(X(idxs(1 : 3, :)), 1); % The actually transmitted 3 symbols
symbols_alt = num2cell(X(idxs(4 : 6, :)), 1); % The alternative 3 symbols

%% 3. Channel settings (Let us try AWGN first)
channel = 'AWGN'; % Channel model, can be specified as AWGN, Rayleigh, Rician, Rayleigh_imp, Rician_imp
if strcmp(channel, 'AWGN') % AWGN channel
    mu_h = [1; 1; 1];
	sigma2_h = [0; 0; 0];
	sigma2_eps = [0; 0; 0];
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

%% 4. Compute PEP for the symbols
N = 64;
xi = 1 / 4;
tic;

PEP_MGF = NaN(1, Q ^ 6); % PEP computed with MGF method

idxs_cell = num2cell(idxs, 1);
interval = floor(Q ^ 6 / 100);
tic;
matlabpool open 4 % My computers has a Core i7-4790 CPU with 4 cores and it takes ~8000 sec. If your computer has more cores it is possible to open more thread to further speed it up
parfor q = 1 : Q ^ 6
    if all(idxs_cell{q}([3,1,2]) ~= idxs_cell{q}([6,5,4]))
        PEP_MGF(q) = get_PEP_symbol(symbols_base{q}, symbols_alt{q}, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi);
    end
end
matlabpool close
toc;


%% 5. Compute PEP for the bits
B = get_n_diff_bits(Nbps);
PEP_MGF_bit = zeros(1, Q ^ 6); % PEP computed with MGF method

for q = 1 : Q ^ 6
    if ~isnan(PEP_MGF(q))
        PEP_MGF_bit(q) = PEP_MGF(q) * B(idxs_cell{q}(1), idxs_cell{q}(2)) / Q;
    end
end

%% 6. Write the result to files using the same format as in http://www.seas.upenn.edu/~hahn/#E
filename = ['q3aptest_', num2str(Q), type_mod, '_', num2str(Eb2N0), 'dB.data'];
fileID = fopen(filename, 'w+');
tic;
fprintf(fileID, '  %18.16e', PEP_MGF_bit);
toc;
fclose(fileID);
