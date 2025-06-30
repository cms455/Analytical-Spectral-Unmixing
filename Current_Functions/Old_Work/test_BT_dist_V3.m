min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 0, 0];  % Change this freely
num_points = 120;
wavelengths = linspace(min_w, max_w, num_points);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
%A = generate_spectrum_curve(num_points, 3, 200,400,2);
k = 6;  % You can change this
delta = 5;
num_iters = 300;
num_columns = size(A, 2);  % Adjust to actual number of columns in A
selection_counts = zeros(1, num_columns);

max_time = 3; % seconds
h = waitbar(0, 'Running iterations...');

% --- Run iterations and accumulate selected indices ---
for i = 1:num_iters
    waitbar(i / num_iters, h, sprintf('Iteration %d of %d', i, num_iters));
    try
        tic;
        [conditioned_indices, submatrix] = bourgain_tzafriri_v2(A, delta);
        elapsed_time = toc;

        if elapsed_time > max_time
            fprintf("Iteration %d skipped (%.2f sec > %.2f sec)\n", i, elapsed_time, max_time);
            continue;
        end

        selection_counts(conditioned_indices) = selection_counts(conditioned_indices) + 1;

    catch ME
        fprintf("Iteration %d failed: %s\n", i, ME.message);
    end
end
close(h);

% --- Find and rank peaks ---
[~, peak_locs] = findpeaks(selection_counts);
peak_frequencies = selection_counts(peak_locs);
[sorted_freqs, sort_idx] = sort(peak_frequencies, 'descend');
sorted_peak_indices = peak_locs(sort_idx);

num_top_peaks = min(10, length(sorted_peak_indices));
top_peak_indices = sorted_peak_indices(1:num_top_peaks);
top_peak_frequencies = sorted_freqs(1:num_top_peaks);

% --- Find valleys (local minima) ---
[~, valley_locs] = findpeaks(-selection_counts);
valley_frequencies = selection_counts(valley_locs);

% --- Always include index 1 as a valley
valley_locs = union([1], valley_locs);  % Always include 1
valley_frequencies = selection_counts(valley_locs);

% --- Check and include last index if it's a valley
last_idx = length(selection_counts);
if selection_counts(end) < selection_counts(end - 1)
    valley_locs = [valley_locs, last_idx];
    valley_frequencies = [valley_frequencies, selection_counts(end)];
end

% --- Sort valleys by ascending frequency
[sorted_valley_freqs, valley_sort_idx] = sort(valley_frequencies, 'ascend');
sorted_valley_indices = valley_locs(valley_sort_idx);

num_valleys = min(10, length(sorted_valley_indices));
top_valley_indices = sorted_valley_indices(1:num_valleys);
top_valley_frequencies = sorted_valley_freqs(1:num_valleys);

% --- Display Peaks ---
fprintf('Top %d Ranked Peaks (Index : Frequency):\n', num_top_peaks);
for i = 1:num_top_peaks
    fprintf('%3d : %d\n', top_peak_indices(i), top_peak_frequencies(i));
end

% --- Display Valleys ---
fprintf('\nTop %d Ranked Valleys (Index : Frequency):\n', num_valleys);
for i = 1:num_valleys
    fprintf('%3d : %d\n', top_valley_indices(i), top_valley_frequencies(i));
end

% --- Plot selection frequency and mark peaks/valleys ---
figure;
bar(selection_counts, 'FaceColor', [0.2 0.6 0.8]);
xlabel('Column Index');
ylabel('Selection Frequency');
title('Column Selection Frequency Across Iterations');
grid on;
hold on;

% Mark peaks with red stars
plot(peak_locs, selection_counts(peak_locs), 'r*', 'MarkerSize', 8);

% Mark valleys with green circles
plot(valley_locs, selection_counts(valley_locs), 'go', 'MarkerSize', 8);

legend('Selection Frequency', 'Detected Peaks', 'Detected Valleys');

% --- Combine peaks and valleys into one search pool ---
search_pool = union(top_peak_indices, top_valley_indices);
search_pool = union(search_pool, [1, num_columns]);
search_pool = unique(search_pool);  % Ensure no duplicates



if length(search_pool) < k
    error('Not enough unique columns in the combined pool to choose %d.', k);
end

combs = nchoosek(search_pool, k);
num_combs = size(combs, 1);

best_norm = Inf;
best_combo = [];

fprintf('\nSearching all %d combinations of %d columns from peaks + valleys...\n', num_combs, k);

for i = 1:num_combs
    cols = combs(i, :);
    subA = A(:, cols);
    pseudo_inv = pinv(subA);
    fro_norm = norm(pseudo_inv, 'fro');

    if fro_norm < best_norm
        best_norm = fro_norm;
        best_combo = cols;
    end
end

% --- Expand best_combo by including ±3 neighbors of each index ---
n_range = 3;
neighborhood = -n_range:n_range;
expanded_set = [];

for i = 1:length(best_combo)
    idx = best_combo(i);
    neighbors = idx + neighborhood;
    % Keep indices within valid bounds
    neighbors = neighbors(neighbors >= 1 & neighbors <= num_columns);
    expanded_set = [expanded_set, neighbors];
end

% Remove duplicates and sort
expanded_set = unique(expanded_set);

% --- Generate all k-combinations from the expanded set ---
if length(expanded_set) < k
    error('Not enough unique columns in the expanded neighborhood to choose %d.', k);
end

fprintf('\nSearching all combinations of %d columns from expanded best_combo neighborhood...\n', k);
neigh_combs = nchoosek(expanded_set, k);
num_neigh_combs = size(neigh_combs, 1);

best_neigh_norm = Inf;
best_neigh_combo = [];

for i = 1:num_neigh_combs
    cols = neigh_combs(i, :);
    subA = A(:, cols);
    pseudo_inv = pinv(subA);
    fro_norm = norm(pseudo_inv, 'fro');

    if fro_norm < best_neigh_norm
        best_neigh_norm = fro_norm;
        best_neigh_combo = cols;
    end
end

% --- Output the best combination from neighborhood expansion ---
fprintf('\nBest combo from neighborhood expansion (min Frobenius norm = %.8f):\n', best_neigh_norm);
disp(best_neigh_combo);

% --- Compare it to previous best ---
neigh_diff = best_neigh_norm - best_norm;
fprintf("Neighborhood DIFFERENCE: %.8f\n", neigh_diff);

if neigh_diff < 0
    disp("Improved over original best_combo ✅");
elseif neigh_diff > 0
    disp("Worse than original best_combo ❌");
else
    disp("No change from original best_combo.");
end


% --- Output the best combination ---
fprintf('\nBest combination of %d columns (min Frobenius norm = %.8f):\n', k, best_norm);
disp(best_combo);

[l_submatrix, l_indices] = luke_algorithm(A', k);
luke_norm = norm(pinv(l_submatrix), 'Fro');

A_norm = normalize_columns(A);
num_iters_BT = 10000;
[~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters_BT);

fprintf('\nMinimum Frobenius norm from BT all-fix: %.8f\n', min_inv_val);
fprintf('Best indices (BT all-fix): ');
fprintf('%d ', min_inv_indices);
fprintf('\n');

fprintf('\nMinimum Frobenius norm from Luke: %.8f\n', luke_norm);
fprintf('Best indices (Luke): ');
fprintf('%d ', l_indices);
fprintf('\n');


min_inv_diff = min_inv_val - best_norm;
fprintf("DIFFERENCE: %d",min_inv_diff);
if min_inv_diff < 0
    disp("\n Bad");
elseif min_inv_diff > 0 
    disp("\n Good");
end

min_inv_luke_diff = luke_norm - best_norm;
fprintf("Luke DIFFERENCE: %d",min_inv_luke_diff);
if min_inv_luke_diff < 0
    disp("\n Bad");
elseif min_inv_luke_diff > 0 
    disp("\n Good");
end

n_min_inv_diff = min_inv_val - best_neigh_norm;
fprintf("Neighbor DIFFERENCE: %d",n_min_inv_diff);
if n_min_inv_diff < 0
    disp("\n Bad");
elseif n_min_inv_diff > 0 
    disp("\n Good");
end

n_min_inv_luke_diff = luke_norm - best_neigh_norm;
fprintf("Neighbor Luke DIFFERENCE: %d",n_min_inv_luke_diff);
if n_min_inv_luke_diff < 0
    disp("\n Bad");
elseif n_min_inv_luke_diff > 0 
    disp("\n Good");
end

