clear all;
close all;
clc;

%% 1. Generate the Gray mapped constellation
Nbps = 3;
type_mod = 'QAM';
pwr = 1;
X = get_constellation(Nbps, type_mod, pwr);
Q = 2 ^ Nbps;

%% 2. Channel settings (Let us try AWGN first)
mu_h = [1; 1; 1];
sigma2_h = [0; 0; 0];
sigma2_eps = [0; 0; 0];
    
Eb2N0 = -2; % Eb/N0 in dB
xi = 1 / 4;
N = 100; % Number of serial expansion
sigma2_v = pwr ./ (Nbps * 10 .^ (Eb2N0 / 10));

map = [1 : Q; [3, 2, 5, 4, 7, 6, 1, 0] + 1];

%% 3. Compute the BER upperbound following the same procedure as the Get_C program
BER_upperbound = get_BER_upper_bound(X, map, mu_h, sigma2_h, sigma2_eps, sigma2_v, N, xi)

%% 4. Compute the BER directly by reading the cost matrix file directly
C = dlmread('q3aptest_8QAM_-2dB.data');
cost_solution = get_BER_upper_bound_from_C(map, C)
