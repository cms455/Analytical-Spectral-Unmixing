load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);
Nt = 4;
max_iter = 100;

Nt = 3;

A_inv = pinv(A);

[U,S,V] = svd(A_inv);

S_vec = diag(S);

S_vec_sqr = S_vec.^2;

S_vec_sqr(S_vec_sqr == 0) = inf;

[sorted_vec, sorted_S_indices] = sort(S_vec_sqr);

selected_S_indices = sorted_S_indices(1:Nt);

A(:,selected_S_indices);




