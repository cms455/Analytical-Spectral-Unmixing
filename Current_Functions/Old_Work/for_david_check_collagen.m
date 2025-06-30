spectrum_matrix = [4345, 1548, 761, 723, 769, 778,693,359,206 ; 442, 586, 816, 856, 1209,1220, 1214, 1128,1024];


A_norm = normalize_columns(spectrum_matrix);
k = 2;
num_iters = 1000;
[~, min_inv_indices, ~, ~, ~, min_inv_val]=bourgain_tzafriri_all_fix_selections(spectrum_matrix, A_norm, k, num_iters);

[l_submatrix, l_indices] = luke_algorithm(spectrum_matrix', k);

disp(min_inv_indices)

disp(l_indices)