min_w = 680;
max_w = 1000;
species_bool = [1, 1, 0, 0, 0];  % Change this freely
num_points = 150;
wavelengths = linspace(min_w, max_w, num_points);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);

delta = 0.8;
num_iters = 1000;
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

        % Count how often each column index is selected
        selection_counts(conditioned_indices) = selection_counts(conditioned_indices) + 1;

    catch ME
        fprintf("Iteration %d failed: %s\n", i, ME.message);
    end
end

close(h);

% --- Find peaks in the selection counts ---
[~, peak_locs] = findpeaks(selection_counts);

% Sort peaks by how often they were selected
[~, sort_idx] = sort(selection_counts(peak_locs), 'descend');
sorted_peak_indices = peak_locs(sort_idx);

% Return top 1/4 of the peak indices
num_top = ceil(0.25 * length(sorted_peak_indices));
top_peak_indices = sorted_peak_indices(1:num_top);

% Display result
disp("Top 1/4 most frequent peak indices:");
disp(top_peak_indices);

% Optional: Plot selection count histogram
figure;
bar(selection_counts, 'FaceColor', [0.2 0.6 0.8]);
xlabel('Column Index');
ylabel('Selection Frequency');
title('Column Selection Frequency Across Iterations');
grid on;

% Optional: Mark peaks
hold on;
plot(peak_locs, selection_counts(peak_locs), 'r*', 'MarkerSize', 8);
legend('Selection Frequency', 'Detected Peaks');

% --- Find peaks in the selection counts ---
[~, peak_locs] = findpeaks(selection_counts);

% Get frequencies at those peaks
peak_frequencies = selection_counts(peak_locs);

% Sort peaks by frequency (descending)
[sorted_freqs, sort_idx] = sort(peak_frequencies, 'descend');
sorted_peak_indices = peak_locs(sort_idx);

num_top = min(10, length(sorted_peak_indices));  % Cap at available number of peaks
top_peak_indices = sorted_peak_indices(1:num_top);
top_peak_frequencies = sorted_freqs(1:num_top);

% Display results as ranked list
fprintf('Top %d Ranked Peaks (Index : Frequency):\n', num_top);
for i = 1:num_top
    fprintf('%3d : %d\n', top_peak_indices(i), top_peak_frequencies(i));
end

% Choose top 10 peaks for brute-force search
pool_size = 15;
top_k_search = min(pool_size, length(sorted_peak_indices)); % just in case < 10
peak_pool = sorted_peak_indices(1:top_k_search);

k = 3;  % number of columns in each combination (change as needed)

combs = nchoosek(peak_pool, k);
num_combs = size(combs, 1);

best_norm = Inf;
best_combo = [];

fprintf('\nSearching all %d combinations of %d columns...\n', num_combs, k);

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

% --- Output the best combination ---
fprintf('\nBest combination of %d columns (min Frobenius norm = %.4f):\n', k, best_norm);
disp(best_combo);

A_norm = normalize_columns(A);
num_iters_BT = 10000;
[~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters_BT);

% Print the result
fprintf('Minimum Frobenius norm: %.4f\n', min_inv_val);
fprintf('Best indices: ');
fprintf('%d ', min_inv_indices);
fprintf('\n');
