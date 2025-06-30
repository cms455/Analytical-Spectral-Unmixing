min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 1, 1];  % Modify freely
num_points = 120;
wavelengths = linspace(min_w, max_w, num_points);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
l_num_points = 200;
A_luke = build_absorption_matrix(min_w, max_w, species_bool, l_num_points);
delta = 14;
num_iters = 1000;
num_columns = size(A, 2);
selection_counts = zeros(1, num_columns);

max_time = 3;  % seconds
h = waitbar(0, 'Running iterations...');

% --- Run Bourgain-Tzafriri and accumulate selection counts ---
for i = 1:num_iters
    waitbar(i / num_iters, h, sprintf('Iteration %d of %d', i, num_iters));
    try
        tic;
        [conditioned_indices, ~] = bourgain_tzafriri_v2(A, delta);
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

% --- Compute envelope of selection frequency signal ---
[env_upper, env_lower] = envelope(selection_counts, 10, 'peak');

% --- Detect peaks on the upper envelope ---
[~, peak_locs] = findpeaks(env_upper);
peak_frequencies = selection_counts(peak_locs);
if selection_counts(1) > selection_counts(2)
    peak_locs = [1, peak_locs];
    peak_frequencies = [selection_counts(1), peak_frequencies];
end
if selection_counts(end) > selection_counts(end - 1)
    peak_locs = [peak_locs, num_columns];
    peak_frequencies = [peak_frequencies, selection_counts(end)];
end
[sorted_freqs, sort_idx] = sort(peak_frequencies, 'descend');
sorted_peak_indices = peak_locs(sort_idx);
num_top_peaks = min(10, length(sorted_peak_indices));
top_peak_indices = sorted_peak_indices(1:num_top_peaks);
top_peak_frequencies = sorted_freqs(1:num_top_peaks);

% --- Detect valleys on the lower envelope ---
[~, valley_locs] = findpeaks(-env_lower);
valley_frequencies = selection_counts(valley_locs);
if selection_counts(1) < selection_counts(2)
    valley_locs = [1, valley_locs];
    valley_frequencies = [selection_counts(1), valley_frequencies];
end
if selection_counts(end) < selection_counts(end - 1)
    valley_locs = [valley_locs, num_columns];
    valley_frequencies = [valley_frequencies, selection_counts(end)];
end
[sorted_valley_freqs, valley_sort_idx] = sort(valley_frequencies, 'ascend');
sorted_valley_indices = valley_locs(valley_sort_idx);
num_valleys = min(10, length(sorted_valley_indices));
top_valley_indices = sorted_valley_indices(1:num_valleys);
top_valley_frequencies = sorted_valley_freqs(1:num_valleys);

% --- Print peaks ---
fprintf('\nTop %d Ranked Peaks (Index : Frequency):\n', num_top_peaks);
for i = 1:num_top_peaks
    fprintf('%3d : %d\n', top_peak_indices(i), top_peak_frequencies(i));
end

% --- Print valleys ---
fprintf('\nTop %d Ranked Valleys (Index : Frequency):\n', num_valleys);
for i = 1:num_valleys
    fprintf('%3d : %d\n', top_valley_indices(i), top_valley_frequencies(i));
end

% --- Plot selection frequency + envelopes + peaks/valleys ---
figure;
bar(selection_counts, 'FaceColor', [0.2 0.6 0.8]);
hold on;
plot(env_upper, '--r', 'LineWidth', 1.2);
plot(env_lower, '--g', 'LineWidth', 1.2);
plot(peak_locs, selection_counts(peak_locs), 'r*', 'MarkerSize', 8);
plot(valley_locs, selection_counts(valley_locs), 'go', 'MarkerSize', 8);
xlabel('Column Index');
ylabel('Selection Frequency');
title('Column Selection Frequency with Envelopes');
legend('Selection Count', 'Upper Envelope', 'Lower Envelope', 'Peaks', 'Valleys');
grid on;

% --- Brute-force search: peaks + valleys, pick best k-columns ---
search_pool = union(top_peak_indices, top_valley_indices);
search_pool = unique(search_pool);
k = 5;

if length(search_pool) < k
    error('Not enough unique indices to form combinations of size %d.', k);
end

combs = nchoosek(search_pool, k);
num_combs = size(combs, 1);

best_norm = Inf;
best_combo = [];

fprintf('\nSearching all %d combinations of %d columns from peaks + valleys...\n', num_combs, k);
for i = 1:num_combs
    cols = combs(i, :);
    subA = A(:, cols);
    fro_norm = norm(pinv(subA), 'fro');
    if fro_norm < best_norm
        best_norm = fro_norm;
        best_combo = cols;
    end
end

fprintf('\nBest combination of %d columns (min Frobenius norm = %.4f):\n', k, best_norm);
disp(best_combo);

% --- Run BT all-fix version for comparison ---

[l_submatrix, l_indices] = luke_algorithm(A_luke', k);
luke_norm = norm(pinv(l_submatrix), 'Fro');

A_norm = normalize_columns(A);
num_iters_BT = 10000;
[~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters_BT);

fprintf('\nMinimum Frobenius norm from BT all-fix: %.4f\n', min_inv_val);
fprintf('Best indices (BT all-fix): ');
fprintf('%d ', min_inv_indices);
fprintf('\n');

fprintf('\nMinimum Frobenius norm from Luke: %.4f\n', luke_norm);
fprintf('Best indices (BT all-fix): ');
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

