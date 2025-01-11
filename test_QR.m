load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);
A = normalize_columns(A);
Nt = 4;
max_iter = 100;
% QR decomposition with column pivoting
[Q, R, P] = qr(A, 'vector');

% P is a permutation vector indicating the column order
disp(P);

