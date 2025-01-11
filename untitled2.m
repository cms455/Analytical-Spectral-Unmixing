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
test = A(:,idx);
test_inv = pinv(test);

idx_2 = [1,137];
test_2 = A(:,idx_2);
test_2_inv = inv(test_2);

