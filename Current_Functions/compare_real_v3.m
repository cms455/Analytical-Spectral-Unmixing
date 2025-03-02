% Define wavelength range
min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 1, 1]; % Define active species
num_points = 40;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
species_counts = [2, 3, 4, 5]; % Incremental species counts
k_items = 8;
num_iters = 100000;
num_repeat = 2; % Number of repetitions for each iteration
colors = lines(length(species_counts));

% Initialize figuresbt_std_vals_holder
figure(1);
hold on;
for i = 1:min(species_counts)-1
    plot(wavelengths, A(i,:), 'LineWidth', 2);
end
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('Initial Spectral Absorption Curves');
legend(arrayfun(@(x) sprintf('Species %d', x), 1:min(species_counts)-1, 'UniformOutput', false));
grid on;
hold off;

figure(2);
hold on;

% Storage for results
bt_mean_vals_holder = zeros(length(species_counts), k_items);
luke_mean_vals_holder = zeros(length(species_counts), k_items);
bt_std_vals_holder = zeros(length(species_counts), k_items);
luke_std_vals_holder = zeros(length(species_counts), k_items);

% Loop through increasing numbers of species
for m = 1:length(species_counts)
    s = species_counts(m);
    
    % Extract the first s species from A
    curve = A(1:s, :); 

    % Plot the added spectral curve
    figure(2);
    hold on;
    for j = 1:size(curve, 1)
        plot(wavelengths, curve(j, :), 'LineWidth', 2, 'DisplayName', sprintf('Species %d', j));
    end
    xlabel('Wavelength (nm)');
    ylabel('Absorption');
    title('Incremental Spectral Absorption Curves');
    legend('-DynamicLegend');
    grid on;
    hold off;

    % Normalize A matrix with the selected species
    A_norm = normalize_columns(A(1:s, :));

    num_elems = k_items - s + 1;

    % Storage for results
    bt_search_val_holder = zeros(num_repeat, num_elems);
    bt_times = zeros(num_repeat, num_elems);
    luke_search_val_holder = zeros(1, num_elems);
    luke_times = zeros(1, num_elems);

    % Initialize waitbar
    total_iterations = num_repeat * num_elems;
    current_iteration = 0;
    wb = waitbar(0, 'Running algorithms, please wait...');

    % Precompute Luke's algorithm results for all k
    for k = 1:num_elems
        tic;
        [l_submatrix, ~] = luke_algorithm(A(1:s, :)', s + k - 1);
        luke_search_val_holder(k) = norm(pinv(l_submatrix), 'Fro');
    end

    % Run Bourgain-Tzafriri algorithm and update waitbar
    for k = 1:num_elems
        for r = 1:num_repeat
            tic;
            [~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A(1:s, :), A_norm, s + k - 1, num_iters);
            bt_search_val_holder(r, k) = min_inv_val;

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
    bt_mean_vals_holder(m, 1:num_elems) = bt_mean_vals;
    bt_std_vals_holder(m, 1:num_elems) = bt_std_vals;
    luke_mean_vals_holder(m, 1:num_elems) = luke_search_val_holder;
end

% Final plot for BT and Luke Algorithm results
figure;
hold on;
for i = 1:length(species_counts)
    s = species_counts(i);
    num_elems = k_items - s + 1;
    errorbar(s:k_items, bt_mean_vals_holder(i, 1:num_elems), bt_std_vals_holder(i, 1:num_elems), ...
        'LineWidth', 2, 'DisplayName', sprintf('BT-ALGO (Species = %d)', s), 'Color', colors(i, :));
    plot(s:k_items, luke_mean_vals_holder(i, 1:num_elems), '--', 'LineWidth', 2, ...
        'DisplayName', sprintf('Luke Algorithm (Species = %d)', s), 'Color', colors(i, :));
end
xlabel('k (Wavelength Selections)');
ylabel('Minimum Inverse Frobenius Norm');
title('Comparison of Minimum Inverse Values for BT and Luke Algorithms');
legend('Location', 'Best');
grid on;
hold off;

% Save all necessary data for recreating the figures
save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

idx = randi(9999); % Unique identifier for each run
file_name = ['multi_plot_data_', num2str(idx),'_'];

% Save comparison data
save(fullfile(save_folder, file_name), ...
    'species_counts', 'k_items', 'wavelengths', ...
    'bt_mean_vals_holder', 'bt_std_vals_holder', ...
    'luke_mean_vals_holder', 'A'); 

disp(['All necessary data for recreating the figures has been saved as: ', file_name]);
