min_w = 680;
max_w =970;
species_bool = [1, 1, 0, 0, 0];
num_points = 150;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);

figure;
hold on;
for i = 1:num_species
    plot(A(i,:), 'LineWidth', 2);
end
legend;
A_norm = normalize_columns(A);
k_items = 2;
num_iters = 100000;

num_repeat = 2; % Number of repetitions for each iteration

% Storage for results
bt_search_val_holder = zeros(num_repeat, k_items - 1);
bt_times = zeros(num_repeat, k_items - 1);
bt_selected_wavelengths = cell(num_repeat, k_items - 1); % Store selected wavelengths

% Luke's algorithm storage
luke_search_val_holder = zeros(1, k_items - 1);
luke_times = zeros(1, k_items - 1);
luke_selected_wavelengths = cell(1, k_items - 1); % Store selected wavelengths

% Initialize waitbar
total_iterations = num_repeat * (k_items - 1);
current_iteration = 0;
wb = waitbar(0, 'Running algorithms, please wait...');

% Precompute Luke's algorithm results for all k
for k = 2:k_items
    tic;
    [l_submatrix, l_indices] = luke_algorithm(A', k); % Extract selected indices
    luke_times(k - 1) = toc;
    luke_search_val_holder(k - 1) = norm(pinv(l_submatrix), 'Fro');
    
    % Store selected wavelengths
    luke_selected_wavelengths{k - 1} = wavelengths(l_indices);
end

% Run Bourgain-Tzafriri algorithm and update waitbar
for k = 2:k_items
    for r = 1:num_repeat
        tic;
        [~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters);
        bt_times(r, k - 1) = toc;
        bt_search_val_holder(r, k - 1) = min_inv_val;

        % Store selected wavelengths
        bt_selected_wavelengths{r, k - 1} = wavelengths(min_inv_indices);

        % Update waitbar
        current_iteration = current_iteration + 1;
        waitbar(current_iteration / total_iterations, wb, sprintf('Running algorithms for k = %d, repetition = %d...', k, r));
    end
end

close(wb);

% Calculate mean and standard deviation
bt_mean_vals = mean(bt_search_val_holder, 1);
bt_std_vals = std(bt_search_val_holder, 0, 1);
luke_mean_vals = luke_search_val_holder;
luke_std_vals = zeros(1, k_items - 1);

bt_mean_times = mean(bt_times, 1);
bt_std_times = std(bt_times, 0, 1);
luke_mean_times = luke_times;
luke_std_times = zeros(1, k_items - 1);

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

% Display selected wavelengths
disp('Selected Wavelengths (BT Algorithm):');
for k = 2:k_items
    for r = 1:num_repeat
        fprintf('k = %d, Repeat = %d: ', k, r);
        disp(bt_selected_wavelengths{r, k - 1});
    end
end

disp('Selected Wavelengths (Luke Algorithm):');
for k = 2:k_items
    fprintf('k = %d: ', k);
    disp(luke_selected_wavelengths{k - 1});
end

figure;
hold on;
for i = 1:num_species
    plot(wavelengths, A(i,:), 'LineWidth', 2);
end
xlabel('Wavelength (nm)', 'FontSize', 14);
ylabel('Absorption', 'FontSize', 14);
title('Absorption Spectra with Selected Wavelengths', 'FontSize', 14);
legend('Species 1', 'Species 2', 'Location', 'Best');
grid on;

colors = lines(2); % Generate distinct colors
for k = 2:k_items
    for r = 1:num_repeat
        % BT Algorithm selections
        x_vals_bt = bt_selected_wavelengths{r, k - 1};
        for x = x_vals_bt
            xline(x, '--', 'Color', colors(1,:), 'LineWidth', 1.5);
        end
    end
    % Luke Algorithm selections
    x_vals_luke = luke_selected_wavelengths{k - 1};
    for x = x_vals_luke
        xline(x, '-', 'Color', colors(2,:), 'LineWidth', 1.5);
    end
end

hold off;

% Define the save path
save_path = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data/Hb_data.mat';

% Save the relevant variables
save(save_path, 'wavelengths', 'A', 'bt_mean_vals', 'bt_std_vals', 'luke_mean_vals', 'luke_std_vals', ...
    'bt_mean_times', 'bt_std_times', 'luke_mean_times', 'luke_std_times', 'bt_selected_wavelengths', 'luke_selected_wavelengths', 'k_items');

disp(['Data saved to ' save_path]);
