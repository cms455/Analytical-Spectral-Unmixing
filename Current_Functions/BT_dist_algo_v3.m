function [best_combo, best_norm] = BT_dist_algo_v3(A, k, num_iters, delta, neigh_range, num_neigh_iters)
    num_columns = size(A, 2); 
    selection_counts = zeros(1, num_columns);
    max_time = 3;

    h = waitbar(0, 'Running BT_dist_algo_v2...');
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

    % Final random search in the neighborhood of best_combo
    neighborhood = [];
    for i = 1:length(best_combo)
        idx = best_combo(i);
        neighbors = idx + (-neigh_range:neigh_range);
        neighbors = neighbors(neighbors >= 1 & neighbors <= num_columns);
        neighborhood = [neighborhood, neighbors];
    end
    neighborhood = unique(neighborhood);

    if length(neighborhood) < k
        warning('Neighborhood too small to perform randomized search. Returning current best.');
        return;
    end

    for i = 1:num_neigh_iters
        rand_indices = randsample(neighborhood, k);
        fro_norm = norm(pinv(A(:, rand_indices)), 'fro');
        if fro_norm < best_norm
            best_norm = fro_norm;
            best_combo = rand_indices;
        end
    end
end
