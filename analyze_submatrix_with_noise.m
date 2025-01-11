function [row_means, row_stds, total_mean, total_std] = analyze_submatrix_with_noise(submatrix, num_noise_distributions)
    % Input:
    % submatrix: A 2D numerical array
    % num_noise_distributions: Number of random noise distributions to apply
    %
    % Output:
    % row_means: Mean for each row aggregated across noise distributions
    % row_stds: Standard deviation for each row aggregated across noise distributions
    % total_mean: Mean of all elements aggregated across noise distributions
    % total_std: Standard deviation of all elements aggregated across noise distributions

    % Get matrix dimensions
    [num_rows, num_cols] = size(submatrix);
    
    % Initialize storage for statistics across runs
    row_means_all_runs = zeros(num_rows, num_noise_distributions);  % Row means for each run
    row_stds_all_runs = zeros(num_rows, num_noise_distributions);   % Row standard deviations for each run
    total_means_all_runs = zeros(1, num_noise_distributions);       % Total mean for each run
    total_stds_all_runs = zeros(1, num_noise_distributions);        % Total standard deviation for each run
    
    % Apply noise and calculate statistics for each noise distribution
    for k = 1:num_noise_distributions
        % Generate random noise vector
        N = randn(num_rows, 1);  % Gaussian noise (one value per row)
        P = ones(num_rows,1);
        N = N+ P;
        % Apply noise to the submatrix
        noisy_matrix = submatrix .* N;
        
        % Calculate row-wise statistics for this run
        row_means_all_runs(:, k) = mean(noisy_matrix);  % Row means
        row_stds_all_runs(:, k) = std(noisy_matrix, 0, 2); % Row standard deviations
        
        % Calculate total statistics for this run
        total_means_all_runs(k) = mean(noisy_matrix(:));   % Total mean
        total_stds_all_runs(k) = std(noisy_matrix(:));     % Total standard deviation
    end
    
    % Aggregate row-wise statistics across all runs
    row_means = mean(row_means_all_runs, 2);  % Mean of row means across all runs
    row_stds = mean(row_stds_all_runs, 2);    % Mean of row stds across all runs
    
    % Aggregate total statistics across all runs
    total_mean = mean(total_means_all_runs);  % Mean of total means across all runs
    total_std = mean(total_stds_all_runs);    % Mean of total stds across all runs
end
