clear all;
close all;
clc;

addpath('./functions/');
%% 1. Generate the Gray mapped constellation
Nbps = 4;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);

%% 2. Generate all the Q ^ 6 index vectors [p, i, j, q, k, l]' (or equavilently [p, a, c, q, b, d]')
Q = 2 ^ Nbps;
order = [4, 1, 5, 2, 6, 3];

%% 3. Channel settings (Let us try AWGN first)
channel = 'Rician'; % Channel model, can be specified as AWGN, Rayleigh, Rician, Rayleigh_imp, Rician_imp
K = 10;
amp = [1, sqrt(2), 2]; % The relative amplitude of the R-D link versus S-D link. Assume that the phase of the Relay and the Source can be perfetly aligned.
Eb2N0 = [-2 : 4]; % Eb/N0 in dB
N = [100, 200, 200, 300, 400, 400, 400]; % Number of serial expansion
n_amp = length(amp);
n_Eb2N0 = length(Eb2N0);
xi = 1 / 4;

test_cases = struct();
i_case = 0;
for i_amp = 1 : length(amp)
    for i_Eb2N0 = 1 : n_Eb2N0
        i_case = i_case + 1;
        test_cases(i_case).type = 'Rician';
        test_cases(i_case).mu_h = sqrt(K / (K + 1)) * [1; 1; amp(i_amp)];
        test_cases(i_case).sigma2_h = 1 / (K + 1) * [1; 1; amp(i_amp) ^ 2];
        test_cases(i_case).sigma2_eps = zeros(3, 1);
        test_cases(i_case).Nbps = Nbps;
        test_cases(i_case).constellation = X;
        test_cases(i_case).Eb2N0 = Eb2N0(i_Eb2N0);
        test_cases(i_case).sigma2_v = 1 / (Nbps * 10 ^ (Eb2N0(i_Eb2N0) / 10));
        test_cases(i_case).N = N(i_Eb2N0);
        test_cases(i_case).xi = xi;
        test_cases(i_case).saved_file = ['q3aptest_', num2str(Q), 'QAM_', test_cases(i_case).type, '_case', num2str(i_case), '.data'];
    end
end
n_case = i_case;

%% 4. Compute PEP for the symbols
tic;
matlabpool open 4 % My computers has a Core i7-4790 CPU with 4 cores and it takes ~8000 sec. If your computer has more cores it is possible to open more thread to further speed it up
for i_case = 1 : n_case
    tic;
    sigma2_v_tmp = test_cases(i_case).sigma2_v;
    N_tmp = test_cases(i_case).N;
    Nbps = test_cases(i_case).Nbps;
    Q = 2 ^ Nbps;
    X = test_cases(i_case).constellation;
    mu_h = test_cases(i_case).mu_h;
    sigma2_h = test_cases(i_case).sigma2_h;
    sigma2_eps = test_cases(i_case).sigma2_eps;
    xi = test_cases(i_case).xi;
    
    B = get_n_diff_bits(Nbps);
    
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

    fileID = fopen(test_cases(i_case).saved_file, 'w+');
    fprintf(fileID, '  %18.16e', PEP_MGF_bit);
    toc;
    fclose(fileID);

    disp(['Case ', num2str(i_case), ' completed']);
end
matlabpool close

