load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);


A_norm = normalize_columns(A);
tic;

[conditioned_indices, min_inv_indices, submatrix_cond, submatrix_inv, min_cond_val, min_inv_val] = bourgain_tzafriri_k_cols_n_iters(A,A_norm, 4,3000);
%[conditioned_indices, min_inv_indices, submatrix_cond, submatrix_inv, min_cond_val, min_inv_val] = bourgain_tzafriri_k_cols_thresh(A,A_norm, 2,1000);
toc;

disp('Conditioned Indices:')
disp(conditioned_indices);

disp('Min Inv Indices:')
disp(min_inv_indices);

A_norm_cond_selected = A_norm(:,conditioned_indices);
A_cond_selected = A(:,conditioned_indices);


A_norm_inv_selected = A_norm(:,min_inv_indices);
A_inv_selected = A(:,min_inv_indices);

fprintf('Normalized Condition Number %d \n',cond(A_norm_cond_selected));
fprintf('Condition Number %d \n',cond(A_cond_selected));
fprintf('Inverse norm: %d \n', norm(pinv(A_cond_selected),'Fro'));


fprintf('Normalized Condition Number %d \n',cond(A_norm_inv_selected));
fprintf('Condition Number %d \n',cond(A_inv_selected));
fprintf('Inverse norm: %d \n', norm(pinv(A_inv_selected),'Fro'));