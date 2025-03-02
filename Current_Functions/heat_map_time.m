% Define range of initial matrix sizes and number of selected wavelengths
min_wavelength = 400; 
max_wavelength = 700;
num_species = 2;

num_wavelengths_list = 50:10:800; % Vary initial matrix size
k_values = 2:8; % Number of selected wavelengths

num_iters = 1000; % Number of iterations for BT Algorithm

% Storage for runtime results
bt_times = zeros(length(num_wavelengths_list), length(k_values));
luke_times = zeros(length(num_wavelengths_list), length(k_values));

total_iterations = length(num_wavelengths_list) * length(k_values);
current_iteration = 0;
wb = waitbar(0, 'Running algorithms, please wait...');


% Loop over different numbers of wavelengths in the initial matrix
for i = 1:length(num_wavelengths_list)
    % Update waitbar progress
    
    num_wavelengths = num_wavelengths_list(i);
    wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);

    % Generate random spectral matrix
    A = rand(num_species, num_wavelengths);
    A_norm = normalize_columns(A);

    % Loop over different numbers of selected wavelengths
    for j = 1:length(k_values)
        current_iteration = current_iteration + 1;
        waitbar(current_iteration / total_iterations, wb, ...
        sprintf('Processing: %d/%d', current_iteration, total_iterations));

        k = k_values(j);
        
        % ---- Bourgain-Tzafriri Algorithm ----
        tic;
        [~, ~, ~, ~, ~, ~] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters);
        bt_times(i, j) = toc;

        % ---- Luke's Algorithm ----
        tic;
        [~, ~] = luke_algorithm(A', k);
        luke_times(i, j) = toc;
    end
end
close(wb);

% ---- Create Heatmap Plots ----
figure;

% Luke Algorithm Heatmap
subplot(1,2,1);
imagesc(num_wavelengths_list, k_values, luke_times'); % Transpose for correct axis orientation
colorbar;
xlabel('Number of Wavelengths in Initial Matrix');
ylabel('Number of Selected Wavelengths (k)');
title('Luke Algorithm Runtime (sec)');
set(gca, 'YDir', 'normal'); % Ensure correct y-axis order

% BT Algorithm Heatmap
subplot(1,2,2);
imagesc(num_wavelengths_list, k_values, bt_times'); % Transpose for correct axis orientation
colorbar;
xlabel('Number of Wavelengths in Initial Matrix');
ylabel('Number of Selected Wavelengths (k)');
title('Bourgain-Tzafriri Algorithm Runtime (sec)');
set(gca, 'YDir', 'normal'); % Ensure correct y-axis order

% Save runtime comparison data
save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
idx = randi(9999);
file_name = ['BT_Luke_runtime_comparison_', num2str(idx), '.mat'];

save(fullfile(save_folder, file_name), ...
    'num_wavelengths_list', 'k_values', 'bt_times', 'luke_times');

disp(['All necessary data saved as: ', file_name]);
