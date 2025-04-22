min_w = 680;
max_w = 800;
species_bool = [1, 1, 1,1,1];
num_points = 20;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);

figure;
hold on;
for i = 1:num_species
    plot(A(i,:),'LineWidth',2);
end
legend;
A_norm = normalize_columns(A);
k_items = 8;
num_iters = 1000;

num_repeat = 2; % Number of repetitions for each iteration

% Storage for results
bt_search_val_holder = zeros(num_repeat, k_items - 1);
bt_times = zeros(num_repeat, k_items - 1);

% Luke's algorithm storage
luke_search_val_holder = zeros(1, k_items - 1);
luke_times = zeros(1, k_items - 1);

% Initialize waitbar
total_iterations = num_repeat * (k_items - 1); % Total iterations for Bourgain-Tzafriri
current_iteration = 0; % Counter for current progress
wb = waitbar(0, 'Running algorithms, please wait...');

% Precompute Luke's algorithm results for all k
for k = 2:k_items
    % Run Luke's algorithm once per k
    tic;
    [l_submatrix, ~] = luke_algorithm(A', k);
    luke_times(k - 1) = toc;
    luke_search_val_holder(k - 1) = norm(pinv(l_submatrix), 'Fro');
end

% Run Bourgain-Tzafriri algorithm and update waitbar
for k = 2:k_items
    for r = 1:num_repeat
        % Bourgain-Tzafriri Algorithm
        tic;
        [min_inv_indices, min_inv_val]= random_search(A,k,num_iters);
        bt_search_val_holder(r, k) = min_inv_val;bt_times(r, k - 1) = toc;
        bt_search_val_holder(r, k - 1) = min_inv_val;

        % Update waitbar
        current_iteration = current_iteration + 1;
        waitbar(current_iteration / total_iterations, wb, sprintf('Running algorithms for k = %d, repetition = %d...', k, r));
    end
end

% Close the waitbar
close(wb);

% Calculate mean and standard deviation
bt_mean_vals = mean(bt_search_val_holder, 1);
bt_std_vals = std(bt_search_val_holder, 0, 1);
luke_mean_vals = luke_search_val_holder; % No random component, so directly use results
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
