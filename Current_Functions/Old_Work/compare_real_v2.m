% Define wavelength range
min_wavelength = 680;
max_wavelength = 800;
num_wavelengths = 20;

wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);
species_counts = [2, 3, 4, 5]; % Incremental species counts
num_repeat = 2; % Number of repetitions for each iteration

% Initial spectrum matrix with a minimal number of species
initial_species_bool = [1, 1]; % First two species
A = build_absorption_matrix(min_wavelength, max_wavelength, initial_species_bool, num_wavelengths);

% Plot initial spectral curves
figure(1);
hold on;
for i = 1:size(A, 1)
    plot(wavelengths, A(i, :), 'LineWidth', 2);
end
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('Initial Spectral Absorption Curves');
legend(arrayfun(@(x) sprintf('Species %d', x), 1:size(A,1), 'UniformOutput', false));
grid on;
hold off;

% Define algorithm parameters
k_items = 8;
num_iters = 1000;
colors = lines(length(species_counts));

% Storage for results
bt_mean_vals_holder = zeros(length(species_counts), k_items);
luke_mean_vals_holder = zeros(length(species_counts), k_items);
bt_std_vals_holder = zeros(length(species_counts), k_items);

% Loop through increasing numbers of species
for m = 1:length(species_counts)
    s = species_counts(m);
    
    % Generate one additional unique spectrum and add it to A
    new_species_bool = zeros(1, length(initial_species_bool) + 1);
    new_species_bool(1:end-1) = initial_species_bool;
    new_species_bool(end) = 1; % Add one more species
    new_curve = build_absorption_matrix(min_wavelength, max_wavelength, new_species_bool, num_wavelengths);
    new_curve = new_curve(end, :); % Grab only the new curve

    A = cat(1, A, new_curve);
    A_norm = normalize_columns(A);

    % Plot the updated spectral curves with the new species added
    figure(2);
    hold on;
    plot(wavelengths, new_curve, 'LineWidth', 2, 'DisplayName', sprintf('Species %d', s));
    xlabel('Wavelength (nm)');
    ylabel('Absorption');
    title('Incremental Spectral Absorption Curves');
    legend('-DynamicLegend');
    grid on;
    hold off;

    % Determine how many elements to include in calculations
    num_elems = k_items - s + 1;

    % Storage for results
    bt_search_val_holder = zeros(num_repeat, num_elems);
    bt_times = zeros(num_repeat, num_elems);

    % Luke's algorithm storage
    luke_search_val_holder = zeros(1, num_elems);
    luke_times = zeros(1, num_elems);

    % Initialize waitbar
    total_iterations = num_repeat * num_elems;
    current_iteration = 0;
    wb = waitbar(0, 'Running algorithms, please wait...');

    % Precompute Luke's algorithm results for all k
    for k = 1:num_elems
        tic;
        [l_submatrix, ~] = luke_algorithm(A', s + k - 1);
        luke_search_val_holder(k) = norm(pinv(l_submatrix), 'Fro');
    end

    % Run Bourgain-Tzafriri algorithm and update waitbar
    for k = 1:num_elems
        for r = 1:num_repeat
            tic;
            [~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, s + k - 1, num_iters);
            bt_search_val_holder(r, k) = min_inv_val;

            % Update waitbar
            current_iteration = current_iteration + 1;
            waitbar(current_iteration / total_iterations, wb, sprintf('Running algorithms for k = %d, repetition = %d...', k, r));
        end
    end

    % Close the waitbar
    close(wb);

    % Store mean and std results
    bt_mean_vals = mean(bt_search_val_holder, 1);
    bt_std_vals = std(bt_search_val_holder, 0, 1);

    bt_mean_vals_holder(m, 1:num_elems) = bt_mean_vals;
    bt_std_vals_holder(m, 1:num_elems) = bt_std_vals;
    luke_mean_vals_holder(m, 1:num_elems) = luke_search_val_holder;

    % Update initial_species_bool for the next iteration
    initial_species_bool = new_species_bool;
end

% Plot final comparison of BT and Luke's algorithm
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
