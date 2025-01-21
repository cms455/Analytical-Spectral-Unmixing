% Generate wavelength range
min_wavelength = 400; % Example value
max_wavelength = 700; % Example value
num_species = 3;
num_wavelengths = 200;
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);

% Generate the spectrum matrix A
A = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, 3);
A_norm = normalize_columns(A);
k_items = 6;
num_iters = 200000;
num_repeat = 4; % Number of repetitions for each iteration

% Storage for results
bt_search_val_holder = zeros(num_repeat, k_items - 1);
luke_search_val_holder = zeros(num_repeat, k_items - 1);
bt_times = zeros(num_repeat, k_items - 1);
luke_times = zeros(num_repeat, k_items - 1);

for k = 2:k_items
    for r = 1:num_repeat
        % Bourgain-Tzafriri Algorithm
        tic;
        [~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters);
        bt_times(r, k - 1) = toc;
        bt_search_val_holder(r, k - 1) = min_inv_val;

        % Luke's Algorithm
        tic;
        [l_submatrix, ~] = luke_algorithm(A', k);
        luke_times(r, k - 1) = toc;
        luke_search_val_holder(r, k - 1) = norm(pinv(l_submatrix), 'Fro');
    end
end

% Calculate mean and standard deviation
bt_mean_vals = mean(bt_search_val_holder, 1);
bt_std_vals = std(bt_search_val_holder, 0, 1);
luke_mean_vals = mean(luke_search_val_holder, 1);
luke_std_vals = std(luke_search_val_holder, 0, 1);

bt_mean_times = mean(bt_times, 1);
bt_std_times = std(bt_times, 0, 1);
luke_mean_times = mean(luke_times, 1);
luke_std_times = std(luke_times, 0, 1);

% Plot minimum inverse values comparison with error bars
figure;
hold on;
errorbar(num_species:k_items, bt_mean_vals(num_species-1:end), bt_std_vals(num_species-1:end), '-o', 'DisplayName', 'Bourgain-Tzafriri', 'LineWidth', 2);
errorbar(num_species:k_items, luke_mean_vals(num_species-1:end), luke_std_vals(num_species-1:end), '-o', 'DisplayName', 'Luke Algorithm', 'LineWidth', 2);
hold off;
xlabel('k (Wavelength Selections)', 'FontSize', 14);
ylabel('Minimum Inverse Frobenius Norm', 'FontSize', 14);
title('Comparison of Minimum Inverse Values for BT and Luke Algorithms', 'FontSize', 14);
legend('Location', 'Best', 'FontSize', 12);
grid on;

% Plot runtime comparison with error bars
figure;
hold on;
errorbar(num_species:k_items, bt_mean_times(num_species-1:end), bt_std_times(num_species-1:end), '-o', 'DisplayName', 'Bourgain-Tzafriri', 'LineWidth', 2);
errorbar(num_species:k_items, luke_mean_times(num_species-1:end), luke_std_times(num_species-1:end), '-o', 'DisplayName', 'Luke Algorithm', 'LineWidth', 2);
hold off;
xlabel('k (Wavelength Selections)', 'FontSize', 14);
ylabel('Runtime (seconds)', 'FontSize', 14);
title('Comparison of Runtime for BT and Luke Algorithms', 'FontSize', 14);
legend('Location', 'Best', 'FontSize', 12);
grid on;
