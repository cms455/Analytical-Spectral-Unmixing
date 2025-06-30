%% Spectral Algorithm Comparison Study
% This script compares different algorithms for optimal wavelength selection
% in spectral unmixing applications by minimizing the inverse Frobenius norm.
%
% Algorithms compared:
% 1. Brute Force (exhaustive search) - baseline
% 2. Random Search (Monte Carlo approach)
% 3. Luke Algorithm (greedy selection)
% 4. NK Algorithm (neighborhood search)

%% Setup Parameters
min_wavelength = 400;
max_wavelength = 700;
num_species = 2;
num_wavelengths = 50;
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);

k_items = 5;  % Maximum number of wavelengths to select
x_idx = 1:num_wavelengths;
num_iters = 10000;  % Iterations for random search
num_repeats = 8;    % Number of different spectra to test

% Preallocate storage arrays
results = struct();
results.brute_force = zeros(num_repeats, k_items - 1);
results.random_search = zeros(num_repeats, k_items - 1);
results.luke_algorithm = zeros(num_repeats, k_items - 1);
results.nk_algorithm = zeros(num_repeats, k_items - 1);

% Runtime storage
runtimes = struct();
runtimes.brute_force = zeros(num_repeats, k_items - 1);
runtimes.random_search = zeros(num_repeats, k_items - 1);
runtimes.luke_algorithm = zeros(num_repeats, k_items - 1);
runtimes.nk_algorithm = zeros(num_repeats, k_items - 1);

% Error tracking (difference from brute force optimal)
errors = struct();
errors.random_search = zeros(num_repeats, k_items - 1);
errors.luke_algorithm = zeros(num_repeats, k_items - 1);
errors.nk_algorithm = zeros(num_repeats, k_items - 1);

%% Main Algorithm Comparison Loop
fprintf('Starting algorithm comparison with %d trials...\n', num_repeats);

for trial = 1:num_repeats
    fprintf('\n--- Trial %d/%d ---\n', trial, num_repeats);
    
    % Generate new random spectrum for this trial
    A = generate_spectrum_curve(num_wavelengths, num_species, ...
                               min_wavelength, max_wavelength, 5);
    A_norm = normalize_columns(A);
    
    % Test each k value (number of wavelengths to select)
    for k = 2:k_items
        fprintf('Testing k = %d wavelengths...\n', k);
        
        % Generate all possible combinations for brute force
        combinations = nchoosek(x_idx, k);
        num_combinations = size(combinations, 1);
        
        %% Brute Force Algorithm (Exhaustive Search)
        fprintf('  Running brute force search (%d combinations)...\n', num_combinations);
        tic;
        brute_vals = zeros(1, num_combinations);
        for i = 1:num_combinations
            idx = combinations(i, :);
            brute_vals(i) = norm(pinv(A(:, idx)), 'fro');
        end
        runtimes.brute_force(trial, k - 1) = toc;
        
        [brute_min_val, ~] = min(brute_vals);
        results.brute_force(trial, k - 1) = brute_min_val;
        
        %% Random Search Algorithm
        fprintf('  Running random search (%d iterations)...\n', num_iters);
        tic;
        [~, random_min_val] = random_search(A, k, num_iters);
        runtimes.random_search(trial, k - 1) = toc;
        results.random_search(trial, k - 1) = random_min_val;
        errors.random_search(trial, k - 1) = abs(brute_min_val - random_min_val);
        
        %% Luke Algorithm
        fprintf('  Running Luke algorithm...\n');
        tic;
        [luke_submatrix, ~] = luke_algorithm(A', k);
        runtimes.luke_algorithm(trial, k - 1) = toc;
        luke_norm = norm(pinv(luke_submatrix), 'fro');
        results.luke_algorithm(trial, k - 1) = luke_norm;
        errors.luke_algorithm(trial, k - 1) = abs(brute_min_val - luke_norm);
        
        %% NK Algorithm
        fprintf('  Running NK algorithm...\n');
        tic;
        [~, nk_val] = nk_column_selector(A, k, 1000);
        runtimes.nk_algorithm(trial, k - 1) = toc;
        nk_norm = abs(nk_val);
        results.nk_algorithm(trial, k - 1) = nk_norm;
        errors.nk_algorithm(trial, k - 1) = abs(brute_min_val - nk_norm);
    end
end

%% Generate Summary Statistics
fprintf('\n=== Computing Summary Statistics ===\n');

% Calculate mean and standard deviation across trials
stats = struct();
algorithm_names = {'brute_force', 'random_search', 'luke_algorithm', 'nk_algorithm'};

for i = 1:length(algorithm_names)
    alg_name = algorithm_names{i};
    stats.(alg_name).mean_result = mean(results.(alg_name), 1);
    stats.(alg_name).std_result = std(results.(alg_name), 1);
    stats.(alg_name).mean_runtime = mean(runtimes.(alg_name), 1);
    stats.(alg_name).std_runtime = std(runtimes.(alg_name), 1);
end

% Calculate error statistics (excluding brute force)
error_algorithms = {'random_search', 'luke_algorithm', 'nk_algorithm'};
for i = 1:length(error_algorithms)
    alg_name = error_algorithms{i};
    stats.(alg_name).mean_error = mean(errors.(alg_name), 1);
    stats.(alg_name).std_error = std(errors.(alg_name), 1);
    stats.(alg_name).max_error = max(errors.(alg_name), [], 1);
end

%% Visualization
create_comparison_plots(results, runtimes, errors, stats, k_items, num_repeats);

%% Save Results
save_results(results, runtimes, errors, stats, num_species, k_items, wavelengths);

fprintf('\nAnalysis complete!\n');

%% Function Definitions

function create_comparison_plots(results, runtimes, errors, stats, k_items, num_repeats)
    %% Plot 1: Algorithm Performance Comparison (Mean ± Std)
    figure('Position', [100, 100, 1200, 400]);
    
    subplot(1, 2, 1);
    k_values = 2:k_items;
    hold on;
    
    % Plot mean results with error bars
    errorbar(k_values, stats.brute_force.mean_result, stats.brute_force.std_result, ...
             '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Brute Force');
    errorbar(k_values, stats.random_search.mean_result, stats.random_search.std_result, ...
             '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Random Search');
    errorbar(k_values, stats.luke_algorithm.mean_result, stats.luke_algorithm.std_result, ...
             '-^', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Luke Algorithm');
    errorbar(k_values, stats.nk_algorithm.mean_result, stats.nk_algorithm.std_result, ...
             '-d', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'NK Algorithm');
    
    xlabel('Number of Selected Wavelengths (k)', 'FontSize', 12);
    ylabel('Inverse Frobenius Norm', 'FontSize', 12);
    title('Algorithm Performance Comparison', 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);
    
    % Plot 2: Runtime Comparison
    subplot(1, 2, 2);
    hold on;
    
    semilogy(k_values, stats.brute_force.mean_runtime, '-o', 'LineWidth', 2, ...
             'MarkerSize', 8, 'DisplayName', 'Brute Force');
    semilogy(k_values, stats.random_search.mean_runtime, '-s', 'LineWidth', 2, ...
             'MarkerSize', 8, 'DisplayName', 'Random Search');
    semilogy(k_values, stats.luke_algorithm.mean_runtime, '-^', 'LineWidth', 2, ...
             'MarkerSize', 8, 'DisplayName', 'Luke Algorithm');
    semilogy(k_values, stats.nk_algorithm.mean_runtime, '-d', 'LineWidth', 2, ...
             'MarkerSize', 8, 'DisplayName', 'NK Algorithm');
    
    xlabel('Number of Selected Wavelengths (k)', 'FontSize', 12);
    ylabel('Runtime (seconds, log scale)', 'FontSize', 12);
    title('Runtime Comparison', 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);
    
    %% Plot 2: Error Analysis
    figure('Position', [150, 150, 800, 600]);
    
    subplot(2, 1, 1);
    hold on;
    
    % Plot individual trial errors (thin lines)
    for trial = 1:num_repeats
        plot(k_values, errors.random_search(trial, :), 'Color', [1 0.7 0.7], 'LineWidth', 0.5);
        plot(k_values, errors.luke_algorithm(trial, :), 'Color', [0.7 0.7 1], 'LineWidth', 0.5);
        plot(k_values, errors.nk_algorithm(trial, :), 'Color', [1 0.8 0.6], 'LineWidth', 0.5);
    end
    
    % Plot mean errors (thick lines)
    plot(k_values, stats.random_search.mean_error, 'r-', 'LineWidth', 3, ...
         'DisplayName', 'Random Search (mean)');
    plot(k_values, stats.luke_algorithm.mean_error, 'b-', 'LineWidth', 3, ...
         'DisplayName', 'Luke Algorithm (mean)');
    plot(k_values, stats.nk_algorithm.mean_error, 'Color', [1 0.6 0.2], 'LineWidth', 3, ...
         'DisplayName', 'NK Algorithm (mean)');
    
    xlabel('Number of Selected Wavelengths (k)', 'FontSize', 12);
    ylabel('Absolute Error from Optimal', 'FontSize', 12);
    title('Error Compared to Brute Force Optimal Solution', 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);
    
    % Subplot 2: Error statistics
    subplot(2, 1, 2);
    bar_data = [stats.random_search.mean_error; stats.luke_algorithm.mean_error; stats.nk_algorithm.mean_error]';
    bar(k_values, bar_data);
    xlabel('Number of Selected Wavelengths (k)', 'FontSize', 12);
    ylabel('Mean Absolute Error', 'FontSize', 12);
    title('Mean Error by Algorithm', 'FontSize', 14);
    legend({'Random Search', 'Luke Algorithm', 'NK Algorithm'}, 'Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);
end

function save_results(results, runtimes, errors, stats, num_species, k_items, wavelengths)
    % Create save directory
    save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end
    
    % Generate unique filename
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    file_name = ['algorithm_comparison_', timestamp, '.mat'];
    
    % Save all data
    save(fullfile(save_folder, file_name), ...
         'results', 'runtimes', 'errors', 'stats', ...
         'num_species', 'k_items', 'wavelengths');
    
    fprintf('Results saved as: %s\n', file_name);
    
    % Also save a summary report
    report_name = ['algorithm_summary_', timestamp, '.txt'];
    save_summary_report(fullfile(save_folder, report_name), stats, k_items);
end

function save_summary_report(filename, stats, k_items)
    fid = fopen(filename, 'w');
    if fid == -1
        warning('Could not create summary report file');
        return;
    end
    
    fprintf(fid, 'ALGORITHM COMPARISON SUMMARY REPORT\n');
    fprintf(fid, '==================================\n\n');
    fprintf(fid, 'Generated: %s\n\n', datestr(now));
    
    k_values = 2:k_items;
    
    % Performance summary
    fprintf(fid, 'MEAN PERFORMANCE (Inverse Frobenius Norm)\n');
    fprintf(fid, 'k\tBrute Force\tRandom Search\tLuke Algorithm\tNK Algorithm\n');
    for i = 1:length(k_values)
        fprintf(fid, '%d\t%.6f\t%.6f\t%.6f\t%.6f\n', k_values(i), ...
                stats.brute_force.mean_result(i), stats.random_search.mean_result(i), ...
                stats.luke_algorithm.mean_result(i), stats.nk_algorithm.mean_result(i));
    end
    
    % Error summary
    fprintf(fid, '\nMEAN ERROR FROM OPTIMAL\n');
    fprintf(fid, 'k\tRandom Search\tLuke Algorithm\tNK Algorithm\n');
    for i = 1:length(k_values)
        fprintf(fid, '%d\t%.6f\t%.6f\t%.6f\n', k_values(i), ...
                stats.random_search.mean_error(i), stats.luke_algorithm.mean_error(i), ...
                stats.nk_algorithm.mean_error(i));
    end
    
    % Runtime summary
    fprintf(fid, '\nMEAN RUNTIME (seconds)\n');
    fprintf(fid, 'k\tBrute Force\tRandom Search\tLuke Algorithm\tNK Algorithm\n');
    for i = 1:length(k_values)
        fprintf(fid, '%d\t%.6f\t%.6f\t%.6f\t%.6f\n', k_values(i), ...
                stats.brute_force.mean_runtime(i), stats.random_search.mean_runtime(i), ...
                stats.luke_algorithm.mean_runtime(i), stats.nk_algorithm.mean_runtime(i));
    end
    
    fclose(fid);
    fprintf('Summary report saved as: %s\n', filename);
end