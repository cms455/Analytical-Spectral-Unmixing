function [best_combo, best_norm] = BT_dist_algo(A, k, num_iters, delta)
    num_columns = size(A, 2); 
    selection_counts = zeros(1, num_columns);
    max_time = 3;

    h = waitbar(0, 'Running iterations...');
    for i = 1:num_iters
        waitbar(i / num_iters, h, sprintf('Iteration %d of %d', i, num_iters));
        try
            tic;
            [conditioned_indices, ~] = bourgain_tzafriri_v2(A, delta);
            elapsed_time = toc;
            if elapsed_time > max_time
                continue;
            end
            selection_counts(conditioned_indices) = selection_counts(conditioned_indices) + 1;
        catch
            continue;
        end
    end
    close(h);

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

    n_range = 3;
    neighborhood = -n_range:n_range;
    expanded_set = [];
    for i = 1:length(best_combo)
        idx = best_combo(i);
        neighbors = idx + neighborhood;
        neighbors = neighbors(neighbors >= 1 & neighbors <= num_columns);
        expanded_set = [expanded_set, neighbors];
    end
    expanded_set = unique(expanded_set);
    if length(expanded_set) < k
        warning('Expanded neighborhood too small, skipping second pass.');
        return;
    end

    neigh_combs = nchoosek(expanded_set, k);
    best_neigh_norm = Inf;
    best_neigh_combo = [];
    for i = 1:size(neigh_combs, 1)
        cols = neigh_combs(i, :);
        fro_norm = norm(pinv(A(:, cols)), 'fro');
        if fro_norm < best_neigh_norm
            best_neigh_norm = fro_norm;
            best_neigh_combo = cols;
        end
    end

    % Return the better of the two
    if best_neigh_norm < best_norm
        best_norm = best_neigh_norm;
        best_combo = best_neigh_combo;
    end
end
