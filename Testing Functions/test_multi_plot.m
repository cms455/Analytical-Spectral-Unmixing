% Parameters
min_wavelength = 400; % Example value
max_wavelength = 700; % Example value
num_wavelengths = 20; % Number of wavelengths
k_items = 8; % Maximum value of k (wavelength selections)
num_trials = 2; % Number of trials for error bars
species_counts = [3, 4, 5]; % Number of species to analyze
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);
num_iter = 100;

% Initialize storage for results
results_mean = zeros(length(species_counts), k_items - 1);
results_std = zeros(length(species_counts), k_items - 1);
luke_values_all = zeros(length(species_counts), k_items - 1);

% Loop over species counts
for s = 1:length(species_counts)
    num_species = species_counts(s);
    bt_values = zeros(num_trials, k_items - 1); % Store results for each trial
    
    % Generate the spectrum matrix A (once for Luke's algorithm)
    A_luke = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, 1);
    
    for trial = 1:num_trials
        % Generate the spectrum matrix A
        A = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, trial);
        A_norm = normalize_columns(A);
        
        for k = 2:k_items
            if k >= num_species % Only calculate when k > num_species
                % Run the BT algorithm and store results
                [~, ~, ~, submatrix_inv, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iter);
                bt_values(trial, k - 1) = min_inv_val;
            end
        end
    end
    
    % Compute mean and standard deviation for BT algorithm
    results_mean(s, :) = mean(bt_values, 1);
    results_std(s, :) = std(bt_values, 0, 1);

    % Run Luke's algorithm (only once since it is deterministic)
    for k = 2:k_items
        if k > num_species % Only calculate when k > num_species
            [l_submatrix, ~] = luke_algorithm(A_luke', k);
            luke_values_all(s, k - 1) = norm(pinv(l_submatrix), 'fro');
        end
    end
end

% Plot the results
figure;
hold on;
set(gca, 'FontSize', 14);

colors = lines(length(species_counts)); % Generate distinct colors

for s = 1:length(species_counts)
    num_species = species_counts(s);
    
    % Plot BT algorithm results
    errorbar(2:k_items, results_mean(s, :), results_std(s, :), '-o', ...
        'LineWidth', 2, 'DisplayName', sprintf('BT: %d Species', num_species), ...
        'Color', colors(s, :));
    
    % Plot Luke algorithm results
    plot(2:k_items, luke_values_all(s, :), '--x', ...
        'LineWidth', 2, 'DisplayName', sprintf('Luke: %d Species', num_species), ...
        'Color', colors(s, :));
end

% Customize the plot
xlabel('k (Wavelength Selections)', 'FontSize', 14);
ylabel('Minimum Inverse Frobenius Norm', 'FontSize', 14);
title('Comparison of BT and Luke Algorithms with Varying Species', 'FontSize', 16);
legend('Location', 'Best', 'FontSize', 12);
grid on;
hold off;
