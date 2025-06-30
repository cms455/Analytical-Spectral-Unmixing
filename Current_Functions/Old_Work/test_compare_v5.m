% Generate wavelength range
min_wavelength = 400; % Example value
max_wavelength = 700; % Example value
num_species = 4;
num_wavelengths = 20;
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);

% Generate the spectrum matrix A
A = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, 3);
A_norm = normalize_columns(A);
k_items = 8;
num_iters = 10000;
figure;
hold on;
for i = 1:num_species
    plot(A(i,:),'LineWidth',2);
end
legend;

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
        [~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters);
        bt_times(r, k - 1) = toc;
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
%title('Comparison of Minimum Inverse Values for BT and Luke Algorithms', 'FontSize', 14);
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


% Save all necessary data for recreating the figures in a single file
save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';

% Create directory if it doesn't exist
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

% Generate a random file identifier to avoid overwriting
idx = randi(9999);
file_name = ['BT_Luke_comparison_data_w_error', num2str(idx), '.mat'];

% Save all relevant variables into a .mat file
save(fullfile(save_folder, file_name), ...
    'min_wavelength', 'max_wavelength', 'num_species', 'num_wavelengths', 'wavelengths', ... % Wavelength-related data
    'A', 'A_norm', ...                                  % Spectrum matrix and normalized spectrum
    'k_items', 'num_iters', 'num_repeat', ...           % Algorithm parameters
    'bt_search_val_holder', 'bt_times', ...             % Bourgain-Tzafriri results
    'luke_search_val_holder', 'luke_times', ...         % Luke's algorithm results
    'bt_mean_vals', 'bt_std_vals', 'luke_mean_vals', 'luke_std_vals', ... % Results for minimum inverse values
    'bt_mean_times', 'bt_std_times', 'luke_mean_times', 'luke_std_times'); % Results for runtime

disp(['All necessary data for recreating the figures has been saved into ', file_name]);
