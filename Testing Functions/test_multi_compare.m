% Generate wavelength range
min_wavelength = 400; % Example value
max_wavelength = 700; % Example value
num_species = 2;
num_wavelengths = 20;
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);

% Generate the spectrum matrix A
A = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, 2);
A_norm = normalize_columns(A);
k_items = 8;
x_idx = 1:length(A);
num_iters = 50000;
num_repeats = 8;

% Storage for results


luke_diff_holder = zeros(num_repeats, k_items-1);
bt_diff_holder = zeros(num_repeats, k_items-1);

for r = 1:num_repeats
A = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, 4);
A_norm = normalize_columns(A);
brute_search_val_holder = zeros(1, k_items - 1);
bt_search_val_holder = zeros(1, k_items - 1);
luke_search_val_holder = zeros(1, k_items - 1);
brute_times = zeros(1, k_items - 1);
bt_times = zeros(1, k_items - 1);
luke_times = zeros(1, k_items - 1);

for k = 2:k_items
    combinations = nchoosek(x_idx, k);

    % Brute Force Algorithm
    disp('BF TIME:')
    tic;
    holder = zeros(1, length(combinations));
    for i = 1:length(combinations)
        idx = combinations(i, :);
        holder(i) = norm(pinv(A(:, idx)), 'Fro');
    end
    brute_times(k - 1) = toc;

    [min_val, min_idx] = min(holder, [], 'all');
    selected_indices = combinations(min_idx, :);
    disp('Selected Indices BF:');
    disp(selected_indices);
    disp('Inverse Val BF:');
    disp(holder(min_idx));
    brute_min_val = holder(min_idx);
    brute_search_val_holder(k - 1) = brute_min_val;
    
    % Bourgain-Tzafriri Algorithm
    disp('BT ALGO')
    tic;
    [conditioned_indices, min_inv_indices, submatrix_cond, submatrix_inv, min_cond_val, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters);
    bt_times(k - 1) = toc;

    disp('Selected Indices BT:')
    disp(min_inv_indices)
    disp('Inverse Val BT:')
    disp(min_inv_val)
    bt_search_val_holder(k - 1) = min_inv_val;

    % Luke's Algorithm
    disp('Luke ALGO')
    tic;
    [l_submatrix, l_indices] = luke_algorithm(A', k);
    luke_times(k - 1) = toc;

    disp('Selected Indices LUKE:')
    disp(l_indices)
    disp('Inverse Val LUKE:')
    min_luke_val = norm(pinv(l_submatrix), 'Fro');
    disp(min_luke_val)
    luke_search_val_holder(k - 1) = min_luke_val;

    luke_diff_holder(r,k-1) = abs(brute_min_val - min_luke_val);
    bt_diff_holder(r,k-1) = abs(brute_min_val - min_inv_val);
end



% Plot minimum inverse values comparison
figure;
set(gca, 'FontSize',14)
hold on;
plot(num_species:k_items, brute_search_val_holder(num_species-1:end), '-o', 'DisplayName', 'Brute Force', 'LineWidth',3);
plot(num_species:k_items, bt_search_val_holder(num_species-1:end), '-o', 'DisplayName', 'Bourgain-Tzafriri','LineWidth',2);
plot(num_species:k_items, luke_search_val_holder(num_species-1:end), '-o', 'DisplayName', 'Luke Algorithm','LineWidth',2);
hold off;
xlabel('k (Wavelength Selections)');
ylabel('Minimum Inverse Frobenius Norm');
title('Comparison of Minimum Inverse Values for Different Algorithms');
legend('Location', 'Best','FontSize',12);


% Plot runtime comparison
figure;
hold on;
set(gca, 'FontSize',14)
plot(num_species:k_items, brute_times(num_species-1:end), '-o', 'DisplayName', 'Brute Force','LineWidth',2);
plot(num_species:k_items, bt_times(num_species-1:end), '-o', 'DisplayName', 'Bourgain-Tzafriri','LineWidth',2);
plot(num_species:k_items, luke_times(num_species-1:end), '-o', 'DisplayName', 'Luke Algorithm','LineWidth',2);
hold off;
xlabel('k (Wavelength Selections)');
ylabel('Runtime (seconds)');
title('Comparison of Runtime for Different Algorithms');
legend('Location', 'Best');
end

figure;
hold on;
plot(luke_diff_holder','b-', 'LineWidth',2);
plot(bt_diff_holder', 'LineWidth',2 );

% Save all necessary data for recreating the figures
save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

idx = randi(9999); % Unique identifier for each run
file_name = ['brute_vs_all_comparison_data_', num2str(idx),'_'];

% Save comparison data
save(fullfile(save_folder, file_name), ...
    'num_species', 'k_items', 'wavelengths', ...
    'brute_search_val_holder', 'bt_search_val_holder', 'luke_search_val_holder', ...
    'brute_times', 'bt_times', 'luke_times', 'luke_diff_holder', 'bt_diff_holder');

disp(['All necessary data for recreating the figures has been saved as: ', file_name]);
