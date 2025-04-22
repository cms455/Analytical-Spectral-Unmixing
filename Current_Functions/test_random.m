load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR_Spectrum.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);

A = randn(2,150);
num_cols = size(A,2);

k = 10;
num_iters = 1000;
holder_indices = zeros(num_iters,k);
holder_min_inv = zeros(1,num_iters);
for i = 1:num_iters
 indices = squeeze(randi(num_cols,1,k));
 while length(unique(indices)) ~= k
     indices = squeeze(randi(num_cols,1,k)); 
 end
 submatrix = A(:,indices);
 norm_inv = norm(pinv(submatrix),'Fro');
 holder_min_inv(i) = norm_inv;
 holder_indices(i,:) = indices;
end

[min_val,min_idx] = min(holder_min_inv);
disp('Minimum Inverse')
disp(min_val)
disp('Indices')
disp(holder_indices(min_idx,:))

A_norm = normalize_columns(A);
[~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters); % Function call
 
disp('Minimum Inverse BT')
disp(min_inv_val)
disp('Indices BT')
disp(min_inv_indices)

[l_submatrix,l_indices] = luke_algorithm(A',k);
luke_min_inv = norm(pinv(l_submatrix),'Fro');

disp('Minimum Inverse Luke')
disp(luke_min_inv);
disp('Indices Luke')
disp(l_indices);