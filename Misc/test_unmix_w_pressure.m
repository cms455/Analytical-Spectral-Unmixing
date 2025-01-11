load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);
Nt = 4;
max_iter = 100;

idx = [1,144];
E_1 = A(:,idx);
E_1_inv = pinv(test);

idx_2 = [1,137];
E_2 = A(:,idx_2);
E_2_inv = inv(test_2);

C = [0.5,0.5];


P_1 = E_1 .* C;
P_2 = E_2 .* C;

N = randn(2,1);

P_1 = P_1 + N;
P_2 = P_2 + N;

noisy_C_1 = E_1_inv .* P_1;
noisy_C_2 = E_2_inv .* P_2;

C_1_n = noisy_C_1 - C_1;
C_2_n = noisy_C_2 - C_2;

c_1_score = (C_1_n) ./ (C_1);
c_2_score = (C_2_n)./(C_2);

fprintf('C1 Score %g \n', c_1_score);
fprintf('C2 Score %g \n', c_2_score);