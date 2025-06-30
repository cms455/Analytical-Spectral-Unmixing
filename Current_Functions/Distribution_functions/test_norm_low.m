% Parameters
min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 1, 1];
num_points = 120;
k = 3;  % number of columns to select
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

% Envelope smoothing
pad_len = 10;
pad_left = fliplr(selection_counts(1:pad_len));
pad_right = fliplr(selection_counts(end-pad_len+1:end));
sel_padded = [pad_left, selection_counts, pad_right];
[env_upper_padded, env_lower_padded] = envelope(sel_padded, pad_len, 'peak');
env_upper = env_upper_padded(pad_len+1:end-pad_len);
env_lower = env_lower_padded(pad_len+1:end-pad_len);
env_diff = abs(env_upper - env_lower);

% Plot envelope difference
figure;
plot(env_diff);
title('Envelope Difference');

% Detect peaks/valleys in envelope difference
[~, peak_locs] = findpeaks(env_diff);
[~, valley_locs] = findpeaks(-env_diff);

% Add endpoints if they're local peaks/valleys
if selection_counts(1) > selection_counts(2), peak_locs = [1, peak_locs]; end
if selection_counts(end) > selection_counts(end-1), peak_locs = [peak_locs, num_cols]; end
if selection_counts(1) < selection_counts(2), valley_locs = [1, valley_locs]; end
if selection_counts(end) < selection_counts(end-1), valley_locs = [valley_locs, num_cols]; end

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
search_pool = unique([top_peaks, top_valleys]);
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

% Optional neighborhood search
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

% Report result
fprintf('\nBest combo (min Frobenius norm = %.4f):\n', best_norm);
disp(best_combo);

[min_inv_indices,min_inv_val] = random_search(A, k, num_trials);
fprintf('Random Search: %.8f', min_inv_val);
disp(min_inv_indices);