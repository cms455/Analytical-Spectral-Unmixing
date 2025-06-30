% Add this to the top of the script
function compare_all_algos()

min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 1, 1];
num_points = 201;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
load_A = load('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/combined_spectra_v2.mat','combined_spectra');
load_A = load_A.combined_spectra;
full_A = load_A(:,2:end)';
%full_A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
species_counts = [2, 3, 4, 5];
k_items = 8;
num_repeat = 1;
num_iters = 100000;
num_nk_iters = 10000;

colors = lines(length(species_counts));
figure(1); clf;
figure(2); clf; hold on;

for i = 1:min(species_counts)-1
    plot(full_A(i, :), 'LineWidth', 2);
end

luke_mean_vals_holder = zeros(length(species_counts), k_items);
luke_std_vals_holder = zeros(length(species_counts), k_items);
random_mean_vals_holder = zeros(length(species_counts), k_items);
random_std_vals_holder = zeros(length(species_counts), k_items);
nk_mean_vals_holder = zeros(length(species_counts), k_items);
nk_std_vals_holder = zeros(length(species_counts), k_items);
random_time_vals_holder = zeros(length(species_counts), k_items);
nk_time_vals_holder = zeros(length(species_counts), k_items);

A = full_A(1, :);

for m = 1:length(species_counts)
    s = species_counts(m);
    curve = full_A(m + 1, :);
    figure(2); hold on;
    plot(curve, 'LineWidth', 2);
    A = cat(1, A, curve);
    A_norm = normalize_columns(A);

    num_elems = k_items - s + 1;
   
    random_search_val_holder = zeros(num_repeat, num_elems);
    nk_search_val_holder = zeros(num_repeat, num_elems);
    
    random_time_holder = zeros(num_repeat, num_elems);
    nk_time_holder = zeros(num_repeat, num_elems);
    luke_search_val_holder = zeros(1, num_elems);

    wb = waitbar(0, 'Running algorithms, please wait...');
    total_iterations = num_repeat * num_elems;
    current_iteration = 0;

    for k = 1:num_elems
        [l_submatrix, ~] = luke_algorithm(A', s + k - 1);
        luke_search_val_holder(k) = norm(pinv(l_submatrix'), 'fro');
    end

    for k = 1:num_elems
        for r = 1:num_repeat
            t1 = tic;
            [min_inv_indices, min_inv_val] = random_search(A, s + k - 1, num_iters);
            random_search_val_holder(r, k) = min_inv_val;
            random_time_holder(r, k) = toc(t1);
            t3 = tic;
            [nk_indices, nk_val] = nk_column_selector(A, s + k - 1, num_nk_iters);
            nk_search_val_holder(r, k) = abs(nk_val);
            nk_time_holder(r, k) = toc(t3);

            current_iteration = current_iteration + 1;
            waitbar(current_iteration / total_iterations, wb, sprintf('k = %d, rep = %d...', k, r));
        end
    end

    close(wb);

    random_mean_vals_holder(m, 1:num_elems) = mean(random_search_val_holder, 1);
   
    random_mean_vals_holder(m, 1:num_elems) = mean(random_search_val_holder, 1);
    random_std_vals_holder(m, 1:num_elems) = std(random_search_val_holder, 0, 1);
    nk_mean_vals_holder(m, 1:num_elems) = mean(nk_search_val_holder, 1);
    nk_std_vals_holder(m, 1:num_elems) = std(nk_search_val_holder, 0, 1);
   
    random_time_vals_holder(m, 1:num_elems) = mean(random_time_holder, 1);
    nk_time_vals_holder(m, 1:num_elems) = mean(nk_time_holder, 1);
    luke_mean_vals_holder(m, 1:num_elems) = luke_search_val_holder;
end

figure;
hold on;
for i = 1:length(species_counts)
    s = species_counts(i);
    num_elems = k_items - s + 1;
    plot(s:k_items, luke_mean_vals_holder(i, 1:num_elems), '--', 'DisplayName', 'Luke', 'LineWidth', 2, 'Color', colors(i, :));
    errorbar(s:k_items, random_mean_vals_holder(i, 1:num_elems), random_std_vals_holder(i, 1:num_elems), ':', 'LineWidth', 2, 'DisplayName', 'BT Dist', 'Color', colors(i, :));
    errorbar(s:k_items, nk_mean_vals_holder(i, 1:num_elems), nk_std_vals_holder(i, 1:num_elems), '-.', 'LineWidth', 2, 'DisplayName', 'NK Model', 'Color', colors(i, :));
end
legend;

time_figure = figure;
hold on;
title('Algorithm Runtime Comparison');
xlabel('Number of Selected Columns');
ylabel('Time (seconds)');
for i = 1:length(species_counts)
    s = species_counts(i);
    num_elems = k_items - s + 1;
   
    plot(s:k_items, random_time_vals_holder(i, 1:num_elems), 's-', 'DisplayName', 'Random', 'Color', colors(i, :));
    plot(s:k_items, nk_time_vals_holder(i, 1:num_elems), 'd-', 'DisplayName', 'NK Model', 'Color', colors(i, :));
end
legend;
% === Bar Chart Figure for Runtime Comparison (Luke, Random, NK) ===
bar_figure = figure;
hold on;

bar_data = [];       % Each row: [Luke, Random, NK] time
bar_labels = {};     % For xtick labels
luke_time_vals_holder = zeros(length(species_counts), k_items);  % Add this if not already declared

for i = 1:length(species_counts)
    s = species_counts(i);
    num_elems = k_items - s + 1;

    for j = 1:num_elems
        k_val = s + j - 1;

        % Estimate Luke time as a single run (you could measure it more precisely)
        t_luke = tic;
        [l_submatrix, ~] = luke_algorithm(A', k_val);
        luke_time = toc(t_luke);
        luke_time_vals_holder(i, j) = luke_time;

        % Collect all algorithm times for this setting
        r_time = random_time_vals_holder(i, j);
        nk_time = nk_time_vals_holder(i, j);
        bar_data(end + 1, :) = [luke_time, r_time, nk_time];

        bar_labels{end + 1} = sprintf('s%d-k%d', s, k_val);
    end
end

bar(bar_data, 'grouped');
legend({'Luke', 'Random Search', 'NK Algorithm'}, 'Location', 'NorthWest');
xticks(1:length(bar_labels));
xticklabels(bar_labels);
xtickangle(45);
ylabel('Mean Runtime (seconds)');
title('Runtime Comparison: Luke vs Random vs NK');
set(gca, 'FontSize', 12);




save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
idx = randi(9999);
file_name = ['multi_plot_data_with_all_', num2str(idx), '.mat'];
save(fullfile(save_folder, file_name), 'species_counts', 'k_items', 'wavelengths', 'luke_mean_vals_holder', 'random_mean_vals_holder', 'random_std_vals_holder', 'nk_mean_vals_holder', 'nk_std_vals_holder', 'random_time_vals_holder', 'nk_time_vals_holder', 'A');
disp(['Data saved: ', file_name]);
end