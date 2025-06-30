min_w = 680;
max_w = 970;
species_bool = [1, 1, 0, 0, 0];
num_points = 290;
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
k_items = 10;
num_iters = 1000;
delta = 5;
num_neigh_reps = 3;

num_repeat = 1; % Number of repetitions for each iteration

% Storage for results
rand_search_val_holder = zeros(num_repeat, k_items - 1);
rand_search_times = zeros(num_repeat, k_items - 1);
rand_selected_wavelengths = cell(num_repeat, k_items - 1);

btdist_search_val_holder = zeros(num_repeat, k_items - 1);
btdist_search_times = zeros(num_repeat, k_items - 1);
btdist_selected_wavelengths = cell(num_repeat, k_items - 1);

luke_search_val_holder = zeros(1, k_items - 1);
luke_times = zeros(1, k_items - 1);
luke_selected_wavelengths = cell(1, k_items - 1);

% Initialize waitbar
total_iterations = num_repeat * (k_items - 1);
current_iteration = 0;
wb = waitbar(0, 'Running algorithms, please wait...');

% Luke's algorithm results for all k
for k = 2:k_items
    tic;
    [l_submatrix, l_indices] = luke_algorithm(A', k);
    luke_times(k - 1) = toc;
    luke_search_val_holder(k - 1) = norm(pinv(l_submatrix), 'Fro');
    luke_selected_wavelengths{k - 1} = wavelengths(l_indices);
end

% Random Search and BT Distributional
for k = 2:k_items
    for r = 1:num_repeat
        % Random Search
        tic;
        [rand_indices, rand_val] = random_search(A, k, num_iters);
        rand_search_times(r, k - 1) = toc;
        rand_search_val_holder(r, k - 1) = rand_val;
        rand_selected_wavelengths{r, k - 1} = wavelengths(rand_indices);

        % BT Dist Algo
        tic;
        [btdist_indices, btdist_val] = BT_dist_algo_v3(A, k, 300, delta, 3,20000);
        btdist_search_times(r, k - 1) = toc;
        btdist_search_val_holder(r, k - 1) = btdist_val;
        btdist_selected_wavelengths{r, k - 1} = wavelengths(btdist_indices);

        % Update waitbar
        current_iteration = current_iteration + 1;
        waitbar(current_iteration / total_iterations, wb, sprintf('Running k = %d, repeat = %d...', k, r));
    end
end
close(wb);

% Stats
rand_mean_vals = mean(rand_search_val_holder, 1);
rand_std_vals = std(rand_search_val_holder, 0, 1);
rand_mean_times = mean(rand_search_times, 1);
rand_std_times = std(rand_search_times, 0, 1);

btdist_mean_vals = mean(btdist_search_val_holder, 1);
btdist_std_vals = std(btdist_search_val_holder, 0, 1);
btdist_mean_times = mean(btdist_search_times, 1);
btdist_std_times = std(btdist_search_times, 0, 1);

luke_mean_vals = luke_search_val_holder;
luke_std_vals = zeros(1, k_items - 1);
luke_mean_times = luke_times;
luke_std_times = zeros(1, k_items - 1);

% Plot minimum inverse norms
figure;
hold on;
errorbar(num_species:k_items, rand_mean_vals(num_species-1:end), rand_std_vals(num_species-1:end), '-o', 'DisplayName', 'Random Search', 'LineWidth', 2);
errorbar(num_species:k_items, btdist_mean_vals(num_species-1:end), btdist_std_vals(num_species-1:end), '-o', 'DisplayName', 'BT Dist', 'LineWidth', 2);
errorbar(num_species:k_items, luke_mean_vals(num_species-1:end), luke_std_vals(num_species-1:end), '-o', 'DisplayName', 'Luke Algorithm', 'LineWidth', 2);
hold off;
xlabel('k (Wavelength Selections)', 'FontSize', 14);
ylabel('Minimum Inverse Frobenius Norm', 'FontSize', 14);
title('Comparison of Minimum Inverse Values', 'FontSize', 14);
legend('Location', 'Best', 'FontSize', 12);
grid on;

% Plot runtimes
figure;
hold on;
errorbar(num_species:k_items, rand_mean_times(num_species-1:end), rand_std_times(num_species-1:end), '-o', 'DisplayName', 'Random Search', 'LineWidth', 2);
errorbar(num_species:k_items, btdist_mean_times(num_species-1:end), btdist_std_times(num_species-1:end), '-o', 'DisplayName', 'BT Dist', 'LineWidth', 2);
errorbar(num_species:k_items, luke_mean_times(num_species-1:end), luke_std_times(num_species-1:end), '-o', 'DisplayName', 'Luke Algorithm', 'LineWidth', 2);
hold off;
xlabel('k (Wavelength Selections)', 'FontSize', 14);
ylabel('Runtime (seconds)', 'FontSize', 14);
title('Runtime Comparison', 'FontSize', 14);
legend('Location', 'Best', 'FontSize', 12);
grid on;

% Plot selected wavelengths on spectrum
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

colors = lines(3);
for k = 2:5
    for r = 1:num_repeat
        for x = rand_selected_wavelengths{r, k - 1}
            xline(x, '--', 'Color', colors(1,:), 'LineWidth', 1.5);
        end
        for x = btdist_selected_wavelengths{r, k - 1}
            xline(x, '-.', 'Color', colors(2,:), 'LineWidth', 1.5);
        end
    end
    for x = luke_selected_wavelengths{k - 1}
        xline(x, '-', 'Color', colors(3,:), 'LineWidth', 1.5);
    end
end
hold off;

% Save results
%save_path = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data/Hb_data.mat';
%save(save_path, 'wavelengths', 'A', 'rand_mean_vals', 'rand_std_vals', 'btdist_mean_vals', 'btdist_std_vals', 'luke_mean_vals', 'luke_std_vals', 'rand_mean_times', 'rand_std_times', 'btdist_mean_times', 'btdist_std_times', 'luke_mean_times', 'luke_std_times', 'rand_selected_wavelengths', 'btdist_selected_wavelengths', 'luke_selected_wavelengths', 'k_items');
%disp(['Data saved to ' save_path]);
