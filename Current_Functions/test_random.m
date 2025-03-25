load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR_Spectrum.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);

k = 4;
num_iters = 1000;
holder_indices = zeros(num_iters,k);
holder_min_inv = zeros(1,num_iters);
for i = 1:num_iters
 indices = squeeze(randi(150,1,k));
 while length(unique(indices)) ~= k
     indices = squeeze(randi(150,1,k)); 
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