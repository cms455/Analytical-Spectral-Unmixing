% Parameters
min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 1, 1];
num_points = 120;
k = 3;  % submatrix width
num_trials = 5000;

% Build matrix
wavelengths = linspace(min_w, max_w, num_points);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
[num_rows, num_cols] = size(A);

% Step 1: Random submatrix sampling by condition number
cond_vals = zeros(num_trials, 1);
index_sets = zeros(num_trials, k);

for i = 1:num_trials
    idx = randperm(num_cols, k);
    subA = A(:, idx);
    cond_vals(i) = cond(subA);
    index_sets(i, :) = idx;
end

% Step 2: Select top 10% best-conditioned submatrices
[~, sorted_idx] = sort(cond_vals);
top_n = round(0.1 * num_trials);
top_indices = index_sets(sorted_idx(1:top_n), :);

% Step 3: Build frequency distribution of selected columns
flat_idx = top_indices(:);
selection_counts = histcounts(flat_idx, 0.5:1:(num_cols + 0.5));

% Step 4: Envelope smoothing
pad_len = 10;
pad_left = fliplr(selection_counts(1:pad_len));
pad_right = fliplr(selection_counts(end - pad_len + 1:end));
sel_padded = [pad_left, selection_counts, pad_right];
[env_upper_padded, env_lower_padded] = envelope(sel_padded, pad_len, 'peak');
env_upper = env_upper_padded(pad_len+1:end-pad_len);
env_lower = env_lower_padded(pad_len+1:end-pad_len);
env_diff = abs(env_upper - env_lower);

% Step 5: Detect peaks and valleys
[~, peak_locs] = findpeaks(env_diff);
[~, valley_locs] = findpeaks(-env_diff);
if selection_counts(1) > selection_counts(2), peak_locs = [1, peak_locs]; end
if selection_counts(end) > selection_counts(end-1), peak_locs = [peak_locs, num_cols]; end
if selection_counts(1) < selection_counts(2), valley_locs = [1, valley_locs]; end
if selection_counts(end) < selection_counts(end-1), valley_locs = [valley_locs, num_cols]; end

[~, pk_sort] = sort(selection_counts(peak_locs), 'descend');
[~, vl_sort] = sort(selection_counts(valley_locs), 'ascend');
top_peaks = peak_locs(pk_sort(1:min(10, end)));
top_valleys = valley_locs(vl_sort(1:min(10, end)));

% Step 6: Plot frequency and envelope
figure;
bar(selection_counts, 'FaceColor', [0.2 0.6 0.8]);
hold on;
plot(env_upper, '--r', 'LineWidth', 1.2);
plot(env_lower, '--g', 'LineWidth', 1.2);
plot(top_peaks, selection_counts(top_peaks), 'r*', 'MarkerSize', 8);
plot(top_valleys, selection_counts(top_valleys), 'go', 'MarkerSize', 8);
xlabel('Column Index'); ylabel('Selection Count');
title('Top Submatrix Index Frequency (by Condition Number)');
legend('Frequency', 'Upper Env', 'Lower Env', 'Peaks', 'Valleys');
grid on;

% Step 7: Brute-force search over peak/valley combos (by Frobenius norm)
search_pool = unique([top_peaks, top_valleys]);
if length(search_pool) < k
    error('Not enough unique indices to form combinations.');
end

combos = nchoosek(search_pool, k);
best_fro = Inf;
best_combo = [];

for i = 1:size(combos, 1)
    idx = combos(i, :);
    subA = A(:, idx);
    fro = norm(pinv(subA), 'fro');
    if fro < best_fro
        best_fro = fro;
        best_combo = idx;
    end
end

% Step 8: Neighborhood search around best Frobenius solution
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
    subA = A(:, idx);
    fro = norm(pinv(subA), 'fro');
    if fro < best_fro
        best_fro = fro;
        best_combo = idx;
    end
end

% Final result
fprintf('\nBest combo (min Frobenius norm = %.8f):\n', best_fro);
disp(best_combo);


[min_inv_indices,min_inv_val] = random_search(A, k, num_trials);
fprintf('Random Search: %.8f', min_inv_val);
disp(min_inv_indices);