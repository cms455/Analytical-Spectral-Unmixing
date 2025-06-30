function pair_freq_smooth = visualize_pairwise_selection(A, k, num_trials)
% visualize_pairwise_selection: Plots pairwise co-selection heatmap
% for top 10% Frobenius-norm submatrices
%
% Inputs:
%   A           - Input matrix (m x n)
%   k           - Number of columns to select
%   num_trials  - Number of random trials (e.g., 10000)

num_cols = size(A, 2); 
norm_vals = zeros(num_trials, 1);
index_sets = zeros(num_trials, k);

% Step 1: Generate random submatrices
for i = 1:num_trials
    idx = randperm(num_cols, k);
    B = A(:, idx);
    norm_vals(i) = norm(pinv(B), 'fro');
    index_sets(i, :) = idx;
end

% Step 2: Select top 10% lowest norm sets
[~, sorted_idx] = sort(norm_vals);
top_n = round(0.1 * num_trials);
top_indices = index_sets(sorted_idx(1:top_n), :);

% Step 3: Initialize and populate pairwise co-selection matrix
pair_freq = zeros(num_cols, num_cols);
for i = 1:top_n
    idx = top_indices(i, :);
    for j = 1:k
        for l = j+1:k
            pair_freq(idx(j), idx(l)) = pair_freq(idx(j), idx(l)) + 1;
            pair_freq(idx(l), idx(j)) = pair_freq(idx(l), idx(j)) + 1;
        end
    end
end

% Step 4: Plot heatmap
figure;
imagesc(pair_freq);
axis square;
colorbar;
xlabel('Column Index'); ylabel('Column Index');
title(sprintf('Pairwise Co-selection Frequency (Top %.0f%%)', 100 * top_n / num_trials));
set(gca, 'XTick', 1:num_cols, 'YTick', 1:num_cols);


% Step 4b: Smooth using 2D smoothing
pair_freq_smooth = smoothdata(pair_freq, 1, 'gaussian', 5); % smooth along rows
pair_freq_smooth = smoothdata(pair_freq_smooth, 2, 'gaussian', 5); % smooth along cols

% Step 4c: Plot smoothed co-selection heatmap
figure;
imagesc(pair_freq_smooth);
axis square;
colorbar;
xlabel('Column Index'); ylabel('Column Index');
title('Smoothed Pairwise Co-selection Frequency');
end

