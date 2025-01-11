load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);
Nt = 4;
max_iter = 100;

k_items =2 ;
x_idx = 1:length(A);
combinations = nchoosek(x_idx,k_items);

cond_val_holder = zeros(1,size(combinations,1));

for i = 1:size(combinations,1)
    idx = combinations(i,:);
    test_matrix = A(:,idx);
    cond_val = cond(test_matrix);
    cond_val_holder(i) = cond_val;
end

[min_val, min_idx] = min(cond_val_holder);


fprintf('Best combination for 2: %d , %d, %d', combinations(min_idx,1),combinations(min_idx,2));