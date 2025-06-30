% Parameters
min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 1, 1];
num_points = 120;
k = 3;  
num_trials = 10000;

neigh_range = 4;
num_neigh_iters = 1000;
% Generate A
wavelengths = linspace(min_w, max_w, num_points);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
[num_rows, num_cols] = size(A);

% Random sampling
norm_vals = zeros(num_trials, 1);
index_sets = zeros(num_trials, k);

for i = 1:num_trials
    idx = randperm(num_cols, k);
    B = A(:, idx);
    norm_vals(i) = norm(pinv(B), 'fro');
    index_sets(i, :) = idx;
end

% Select top 10% by lowest Frobenius norm
[~, sorted_idx] = sort(norm_vals);
top_n = round(0.1 * num_trials);
top_indices = index_sets(sorted_idx(1:top_n), :);

% Build frequency distribution
flat_idx = top_indices(:);
selection_counts = histcounts(flat_idx, 0.5:1:(num_cols + 0.5));

% Smooth the selection count distribution usi
x = (1:num_cols)';
y = selection_counts(:);

% Choose smoothing method
smooth_span = 0.15; 
y_smooth = smoothdata(y);

% Plot residuals from smoothing
figure;
plot(y_smooth, 'k', 'LineWidth', 1.5);
title('Smoothed Selection Frequency');
xlabel('Column Index'); ylabel('Frequency');

env_diff = y_smooth;

[~, peak_locs, peak_w, peak_prom] = findpeaks(env_diff);  
[~, valley_locs, valley_w, valley_prom] = findpeaks(-env_diff);  
peak_locs = peak_locs';
valley_locs = valley_locs';


% Sort top N peaks/valleys by frequency
[~, pk_sort] = sort(selection_counts(peak_locs), 'descend');
[~, vl_sort] = sort(selection_counts(valley_locs), 'ascend');
top_peaks = peak_locs(pk_sort(1:min(10, end)));
top_valleys = valley_locs(vl_sort(1:min(10, end)));

% Plot frequency and envelopes
figure;
bar(selection_counts, 'FaceColor', [0.2 0.6 0.8]);
hold on;
plot(env_upper, '--r', 'LineWidth', 1.2);
plot(env_lower, '--g', 'LineWidth', 1.2);
plot(top_peaks, selection_counts(top_peaks), 'r*', 'MarkerSize', 8);
plot(top_valleys, selection_counts(top_valleys), 'go', 'MarkerSize', 8);
xlabel('Column Index'); ylabel('Selection Count');
title('Top Submatrix Index Frequency and Envelope Peaks/Valleys');
legend('Frequency', 'Upper Env', 'Lower Env', 'Peaks', 'Valleys');
grid on;

% Combine peak/valley indices and search for best combo
search_pool = unique([top_peaks, top_valleys,1,num_cols]);
if length(search_pool) < k
    error('Not enough indices to form combinations.');
end

combos = nchoosek(search_pool, k);
best_norm = Inf;
best_combo = [];
for i = 1:size(combos, 1)
    idx = combos(i, :);
    B = A(:, idx);
    nrm = norm(pinv(B), 'fro');
    if nrm < best_norm
        best_norm = nrm;
        best_combo = idx;
    end
end

no_neigh_combo = best_combo;

% Neighbor Search
neigh_range = 4;
neigh_iters = 10000;
neighborhood = [];
for i = 1:k
    n = best_combo(i) + (-neigh_range:neigh_range);
    neighborhood = [neighborhood, n(n >= 1 & n <= num_cols)];
end
neighborhood = unique(neighborhood);

for i = 1:neigh_iters
    idx = randsample(neighborhood, k);
    nrm = norm(pinv(A(:, idx)), 'fro');
    if nrm < best_norm
        best_norm = nrm;
        best_combo = idx;
    end
end

% Results
fprintf('\nBest combo (min Frobenius norm = %.8f):\n', best_norm);
disp(best_combo);

rand_num_trials = 100000;
[min_inv_indices,min_inv_val] = random_search(A, k, rand_num_trials);
fprintf('Random Search: %.8f', min_inv_val);
disp(min_inv_indices);


no_neigh_norm = norm(pinv(A(:, no_neigh_combo)),'fro');
fprintf('\nBest combo no NEIGHBOR (min Frobenius norm = %.8f):\n', no_neigh_norm);
disp(no_neigh_combo);

[l_submatrix, l_indices] = luke_algorithm(A', k);
l_norm = norm(pinv(l_submatrix'), 'fro');
fprintf('Luke Search (min Frobenius norm = %.8f):\n', l_norm);
disp(l_indices);
