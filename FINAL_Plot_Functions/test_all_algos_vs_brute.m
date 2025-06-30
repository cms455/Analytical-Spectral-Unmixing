%% Setup Parameters
min_wavelength = 400;
max_wavelength = 700;
num_species = 2;
num_wavelengths = 50;
wavelengths = linspace(min_wavelength, max_wavelength, num_wavelengths);

k_items = 5;
x_idx = 1:num_wavelengths;
num_iters = 2000;
num_repeats = 3;
delta = 10;
num_neigh_reps = 4;

% Preallocate storage for absolute difference to brute force (per trial)
luke_diff_holder = zeros(num_repeats, k_items - 1);
bt_diff_holder = zeros(num_repeats, k_items - 1);
btdist_diff_holder = zeros(num_repeats, k_items - 1);
nk_diff_holder = zeros(num_repeats, k_items - 1);

%% Main Repetition Loop
for r = 1:num_repeats
    % Generate new random spectrum for this trial
    A = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, 5);
    A_norm = normalize_columns(A);

    % Storage for each algorithm this round
    brute_search_val_holder = zeros(1, k_items - 1);
    bt_search_val_holder = zeros(1, k_items - 1);
    btdist_search_val_holder = zeros(1, k_items - 1);
    luke_search_val_holder = zeros(1, k_items - 1);
    nk_search_val_holder = zeros(1, k_items - 1);

    brute_times = zeros(1, k_items - 1);
    bt_times = zeros(1, k_items - 1);
    btdist_times = zeros(1, k_items - 1);
    luke_times = zeros(1, k_items - 1);
    nk_times = zeros(1, k_items - 1);

    for k = 2:k_items
        combinations = nchoosek(x_idx, k);

        %% Brute Force Algorithm
        disp('Brute Force Algorithm')
        tic;
        brute_vals = zeros(1, size(combinations, 1));
        for i = 1:size(combinations, 1)
            idx = combinations(i, :);
            brute_vals(i) = norm(pinv(A(:, idx)), 'Fro');
        end
        brute_times(k - 1) = toc;

        [brute_min_val, min_idx] = min(brute_vals);
        brute_search_val_holder(k - 1) = brute_min_val;

        %% Random Search Algorithm
        disp('Random Search Algorithm')
        tic;
        [~, min_inv_val] = random_search(A, k, num_iters);
        bt_times(k - 1) = toc;
        bt_search_val_holder(k - 1) = min_inv_val;

        %% BT Distributional Search
        disp('BT Distributional Algorithm')
        tic;
        [~, btdist_norm] = og_dist_v1(A, k, 4000, k, 10000);
        btdist_times(k - 1) = toc;
        btdist_search_val_holder(k - 1) = btdist_norm;

        %% Luke Algorithm
        disp('Luke Algorithm')
        tic;
        [l_submatrix, ~] = luke_algorithm(A', k);
        luke_times(k - 1) = toc;
        luke_norm = norm(pinv(l_submatrix), 'Fro');
        luke_search_val_holder(k - 1) = luke_norm;

        %% NK Algorithm
        disp('NK Algorithm')
        tic;
        [~, nk_val] = nk_column_selector(A, k, 1000);
        nk_times(k - 1) = toc;
        nk_search_val_holder(k - 1) = abs(nk_val);

        %% Record Absolute Differences to Brute Force
        luke_diff_holder(r, k - 1) = abs(brute_min_val - luke_norm);
        bt_diff_holder(r, k - 1) = abs(brute_min_val - min_inv_val);
        btdist_diff_holder(r, k - 1) = abs(brute_min_val - btdist_norm);
        nk_diff_holder(r, k - 1) = abs(brute_min_val - nk_search_val_holder(k - 1));
    end

    %% Plot: Inverse Norm Comparison (One Trial)
    figure;
    hold on;
    set(gca, 'FontSize', 14)
    plot(2:k_items, brute_search_val_holder, '-o', 'DisplayName', 'Brute Force', 'LineWidth', 3);
    plot(2:k_items, bt_search_val_holder, '-o', 'DisplayName', 'Random Search', 'LineWidth', 2);
    plot(2:k_items, btdist_search_val_holder, '-o', 'DisplayName', 'BT Dist', 'LineWidth', 2);
    plot(2:k_items, luke_search_val_holder, '-o', 'DisplayName', 'Luke Algorithm', 'LineWidth', 2);
    plot(2:k_items, nk_search_val_holder, '-o', 'DisplayName', 'NK Algorithm', 'LineWidth', 2);
    hold off;
    xlabel('k (Wavelength Selections)');
    ylabel('Minimum Inverse Frobenius Norm');
    title(sprintf('Trial %d: Minimum Inverse Norm Comparison', r));
    legend('Location', 'Best');

    %% Plot: Runtime Comparison (One Trial)
    figure;
    hold on;
    set(gca, 'FontSize', 14)
    plot(2:k_items, brute_times, '-o', 'DisplayName', 'Brute Force', 'LineWidth', 2);
    plot(2:k_items, bt_times, '-o', 'DisplayName', 'Random Search', 'LineWidth', 2);
    plot(2:k_items, btdist_times, '-o', 'DisplayName', 'BT Dist', 'LineWidth', 2);
    plot(2:k_items, luke_times, '-o', 'DisplayName', 'Luke Algorithm', 'LineWidth', 2);
    plot(2:k_items, nk_times, '-o', 'DisplayName', 'NK Algorithm', 'LineWidth', 2);
    hold off;
    xlabel('k (Wavelength Selections)');
    ylabel('Runtime (seconds)');
    title(sprintf('Trial %d: Runtime Comparison', r));
    legend('Location', 'Best');
end

%% Plot: Error Compared to Brute (All Trials + Mean)
figure;
hold on;

% Plot all repetitions
for r = 1:num_repeats
    plot(2:k_items, luke_diff_holder(r, :), 'Color', [0.5 0.5 1], 'LineWidth', 1);
    plot(2:k_items, bt_diff_holder(r, :), 'Color', [1 0.5 0.5], 'LineWidth', 1);
    plot(2:k_items, btdist_diff_holder(r, :), 'Color', [0.5 1 0.5], 'LineWidth', 1);
    plot(2:k_items, nk_diff_holder(r, :), 'Color', [1 0.6 0.2], 'LineWidth', 1);
end

% Plot mean lines across all repetitions
mean_luke_diff = mean(luke_diff_holder, 1);
mean_bt_diff = mean(bt_diff_holder, 1);
mean_btdist_diff = mean(btdist_diff_holder, 1);
mean_nk_diff = mean(nk_diff_holder, 1);
plot(2:k_items, mean_luke_diff, 'b-', 'LineWidth', 3, 'DisplayName', 'Luke (mean)');
plot(2:k_items, mean_bt_diff, 'r-', 'LineWidth', 3, 'DisplayName', 'Random Search (mean)');
plot(2:k_items, mean_btdist_diff, 'g-', 'LineWidth', 3, 'DisplayName', 'BT Dist (mean)');
plot(2:k_items, mean_nk_diff, 'Color', [1 0.6 0.2], 'LineWidth', 3, 'DisplayName', 'NK (mean)');

legend show;
xlabel('k (Wavelength Selections)');
ylabel('Absolute Difference from Brute Force');
title('Error Compared to Brute Force Across Trials');
set(gca, 'FontSize', 14);

%% Save Data
save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

idx = randi(9999);
file_name = ['brute_vs_all_comparison_data_', num2str(idx), '_'];
save(fullfile(save_folder, file_name), ...
    'num_species', 'k_items', 'wavelengths', ...
    'brute_search_val_holder', 'bt_search_val_holder', 'btdist_search_val_holder', 'luke_search_val_holder', 'nk_search_val_holder', ...
    'brute_times', 'bt_times', 'btdist_times', 'luke_times', 'nk_times', ...
    'luke_diff_holder', 'bt_diff_holder', 'btdist_diff_holder', 'nk_diff_holder');
disp(['All data saved as: ', file_name]);
