% Generate wavelength range
min_wavelength = 400; % Example value
max_wavelength = 700; % Example value
num_species = 2;
num_wavelengths = 100;
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);

% Generate the spectrum matrix A
A = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, 7);
A_norm = normalize_columns(A);
k_items = 6;
num_iters = 1000;
delta = 8;

% Storage for results
bt_search_val_holder = zeros(1, k_items - 1);
luke_search_val_holder = zeros(1, k_items - 1);
btdist_search_val_holder = zeros(1, k_items - 1);

bt_times = zeros(1, k_items - 1);
luke_times = zeros(1, k_items - 1);
btdist_times = zeros(1, k_items - 1);

for k = num_species:k_items
    % --- BT Fixed ---
    disp('BT ALGO')
    tic;
    [~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, k, num_iters);
    bt_times(k - 1) = toc;
    bt_search_val_holder(k - 1) = min_inv_val;

    % --- Luke ---
    disp('Luke ALGO')
    tic;
    [l_submatrix, l_indices] = luke_algorithm(A', k);
    luke_times(k - 1) = toc;
    min_luke_val = norm(pinv(l_submatrix), 'fro');
    luke_search_val_holder(k - 1) = min_luke_val;

    % --- BT Distribution ---
    disp('BT Dist ALGO')
    tic;
    [btdist_indices, btdist_norm] = BT_dist_algo_v3(A, k, 300, delta,3,20000);
    btdist_times(k - 1) = toc;
    btdist_search_val_holder(k - 1) = btdist_norm;

    % Plot all rows of A with selected indices
    figure;
    hold on;
    for i = 1:num_species
        plot(wavelengths, A(i, :), 'LineWidth', 2);
    end
    for i = 1:k
        xline(wavelengths(min_inv_indices(i)), 'r--', 'LineWidth', 1.5);
        xline(wavelengths(l_indices(i)), 'b--', 'LineWidth', 1.5);
        xline(wavelengths(btdist_indices(i)), 'g--', 'LineWidth', 1.5);
    end
    hold off;
    xlabel('Wavelength (nm)');
    ylabel('Absorption');
    title(sprintf('k = %d: Absorption Spectra with Selected Indices Highlighted', k));
    legend(arrayfun(@(x) sprintf('Species %d', x), 1:num_species, 'UniformOutput', false), 'Location', 'Best');
end

% --- Inverse Norm Comparison Plot ---
figure; hold on; set(gca,'FontSize',14)
plot(num_species:k_items, bt_search_val_holder(num_species-1:end), '-o', 'DisplayName', 'BT Fix', 'LineWidth',2);
plot(num_species:k_items, luke_search_val_holder(num_species-1:end), '-o', 'DisplayName', 'Luke', 'LineWidth',2);
plot(num_species:k_items, btdist_search_val_holder(num_species-1:end), '-o', 'DisplayName', 'BT Dist', 'LineWidth',2);
xlabel('k (Wavelength Selections)');
ylabel('Minimum Inverse Frobenius Norm');
title('Comparison of Minimum Inverse Values');
legend('Location', 'Best');

% --- Runtime Comparison Plot ---
figure; hold on; set(gca,'FontSize',14)
plot(num_species:k_items, bt_times(num_species-1:end), '-o', 'DisplayName', 'BT Fix', 'LineWidth',2);
plot(num_species:k_items, luke_times(num_species-1:end), '-o', 'DisplayName', 'Luke', 'LineWidth',2);
plot(num_species:k_items, btdist_times(num_species-1:end), '-o', 'DisplayName', 'BT Dist', 'LineWidth',2);
xlabel('k (Wavelength Selections)');
ylabel('Runtime (seconds)');
title('Comparison of Runtime');
legend('Location', 'Best');

% Save all data
save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';
if ~exist(save_folder, 'dir'); mkdir(save_folder); end
idx = randi(9999);
file_name = ['BT_Luke_BTDist_comparison_data_', num2str(idx),'.mat'];
save(fullfile(save_folder, file_name), ...
    'wavelengths', 'A', 'A_norm', 'num_species', 'k_items', ...
    'bt_search_val_holder', 'luke_search_val_holder', 'btdist_search_val_holder', ...
    'bt_times', 'luke_times', 'btdist_times');
disp(['Saved data as ', file_name]);
