% Parameters
min_wavelength = 400; % Example value
max_wavelength = 700; % Example value
num_wavelengths = 30;
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);
species_counts = [3, 4, 5];
num_repeat = 10; % Number of repetitions for BT algorithm
k_items = 8; % Maximum value of k
num_iters = 1000; % Iterations for BT algorithm

% Define color scheme
colors = lines(length(species_counts)); % Generate distinct colors for species counts

% Create a figure for the final comparison
figure;
hold on;
set(gca, 'FontSize', 14);

for k = 1:length(species_counts)
    % Species-specific parameters
    s = species_counts(k);

    % Generate spectrum matrix A
    A = generate_spectrum_curve(num_wavelengths, s, min_wavelength, max_wavelength, k);
    A_norm = normalize_columns(A);

    % Number of k values to consider
    num_elems = k_items - s + 1;

    % Storage for results
    bt_search_val_holder = zeros(num_repeat, num_elems);
    luke_search_val_holder = zeros(1, num_elems);

    % Precompute Luke's algorithm results
    for j = 1:num_elems
        k_val = s + j - 1;
        [l_submatrix, ~] = luke_algorithm(A', k_val);
        luke_search_val_holder(j) = norm(pinv(l_submatrix), 'fro');
    end

    % Run BT algorithm for multiple repetitions
    for j = 1:num_elems
        k_val = s + j - 1;

        for r = 1:num_repeat
            % Bourgain-Tzafriri Algorithm
            [~, ~, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k_val, num_iters);
            bt_search_val_holder(r, j) = min_inv_val;
        end
    end

    % Calculate mean and standard deviation for BT algorithm
    bt_mean_vals = mean(bt_search_val_holder, 1);
    bt_std_vals = std(bt_search_val_holder, 0, 1);

    % Plot BT algorithm results with error bars
    errorbar(s:k_items, bt_mean_vals, bt_std_vals, '-o', ...
        'LineWidth', 2, 'Color', colors(k, :), ...
        'DisplayName', sprintf('BT Algorithm (%d Species)', s));

    % Plot Luke algorithm results
    plot(s:k_items, luke_search_val_holder, '--x', ...
        'LineWidth', 2, 'Color', colors(k, :), ...
        'DisplayName', sprintf('Luke Algorithm (%d Species)', s));
end

% Customize the plot
xlabel('k (Wavelength Selections)', 'FontSize', 14);
ylabel('Minimum Inverse Frobenius Norm', 'FontSize', 14);
title('Comparison of BT and Luke Algorithms for Different Species', 'FontSize', 16);
legend('Location', 'Best', 'FontSize', 12);
grid on;
hold off;
