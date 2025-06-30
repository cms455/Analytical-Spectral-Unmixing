

% --- Main script below this line ---
min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 0, 0];
num_points = 120;
wavelengths = linspace(min_w, max_w, num_points);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
k = 6;
delta = 5;
num_iters = 300;
num_columns = size(A, 2);
selection_counts = zeros(1, num_columns);

max_time = 3;
h = waitbar(0, 'Running iterations...');
for i = 1:num_iters
    waitbar(i / num_iters, h, sprintf('Iteration %d of %d', i, num_iters));
    try
        tic;
        [conditioned_indices, ~] = bourgain_tzafriri_v2(A, delta);
        if toc > max_time
            continue;
        end
        selection_counts(conditioned_indices) = selection_counts(conditioned_indices) + 1;
    catch
        continue;
    end
end
close(h);

[~, peak_locs] = findpeaks(selection_counts);
[~, valley_locs] = findpeaks(-selection_counts);
valley_locs = union([1], valley_locs);
if selection_counts(end) < selection_counts(end - 1)
    valley_locs = [valley_locs, num_columns];
end

[~, peak_locs] = findpeaks(selection_counts);
peak_frequencies = selection_counts(peak_locs);
[~, sort_idx] = sort(peak_frequencies, 'descend');
sorted_peak_indices = peak_locs(sort_idx);
num_top_peaks = min(10, length(sorted_peak_indices));
top_peak_indices = sorted_peak_indices(1:num_top_peaks);

[~, valley_locs] = findpeaks(-selection_counts);
valley_locs = union([1], valley_locs);
if selection_counts(end) < selection_counts(end - 1)
    valley_locs = [valley_locs, num_columns];
end
valley_frequencies = selection_counts(valley_locs);
[~, valley_sort_idx] = sort(valley_frequencies, 'ascend');
sorted_valley_indices = valley_locs(valley_sort_idx);
num_valleys = min(10, length(sorted_valley_indices));
top_valley_indices = sorted_valley_indices(1:num_valleys);

search_pool = unique([top_peak_indices, top_valley_indices, 1, num_columns]);
if length(search_pool) < k
    error('Not enough unique columns in the combined pool to choose %d.', k);
end


combs = nchoosek(search_pool, k);
best_norm = Inf;
best_combo = [];
for i = 1:size(combs, 1)
    cols = combs(i, :);
    fro_norm = norm(pinv(A(:, cols)), 'fro');
    if fro_norm < best_norm
        best_norm = fro_norm;
        best_combo = cols;
    end
end

num_neigh_reps = 4;  % how many times to refine
range = 2;  % neighborhood range
[best_neigh_combo, best_neigh_norm] = repeated_neighborhood_search(A, k, best_combo, range, num_neigh_reps);

fprintf('\nBest combo (norm = %.8f):\n', best_norm);
disp(best_combo);
fprintf('\nBest from neighborhood search (norm = %.8f):\n', best_neigh_norm);
disp(best_neigh_combo);

[l_submatrix, l_indices] = luke_algorithm(A', k);
luke_norm = norm(pinv(l_submatrix), 'fro');

A_norm = normalize_columns(A);
[~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, 10000);


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



% --- Helper function for repeated neighborhood search ---
function [best_combo, best_norm] = repeated_neighborhood_search(A, k, initial_combo, range, reps)
    num_columns = size(A, 2);
    current_combo = initial_combo;
    best_norm = Inf;
    best_combo = [];

    for r = 1:reps
        neighborhood = [];
        for i = 1:length(current_combo)
            idx = current_combo(i);
            neighbors = idx + (-range:range);
            neighbors = neighbors(neighbors >= 1 & neighbors <= num_columns);
            neighborhood = [neighborhood, neighbors];
        end
        neighborhood = unique(neighborhood);

        if length(neighborhood) < k
            warning('Neighborhood too small for k=%d', k);
            break;
        end

        combs = nchoosek(neighborhood, k);
        for i = 1:size(combs, 1)
            cols = combs(i, :);
            subA = A(:, cols);
            fro_norm = norm(pinv(subA), 'fro');
            if fro_norm < best_norm
                best_norm = fro_norm;
                best_combo = cols;
            end
        end
        current_combo = best_combo;
    end
end