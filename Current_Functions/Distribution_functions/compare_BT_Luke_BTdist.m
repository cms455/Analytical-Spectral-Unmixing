% Add this to the top of the script
function compare_BT_Luke_BTdist()

min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 1, 1];
num_points = 120;
%deltas = [4,6,10,14];
deltas = [100,100,200,200];
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
%full_A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
load_A = load('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/combined_spectra_v2.mat','combined_spectra');
load_A = load_A.combined_spectra;
full_A = load_A(:,2:end)';
species_counts = [2, 3, 4, 5];
k_items = 6;
num_repeat = 1;
num_iters = 10000;

%full_A = generate_spectrum_curve(num_points, num_species, min_w, max_w, 3);

colors = lines(length(species_counts));
figure(1); clf;
figure(2); clf; hold on;

for i = 1:min(species_counts)-1
    plot(full_A(i, :), 'LineWidth', 2);
end

bt_mean_vals_holder = zeros(length(species_counts), k_items);
luke_mean_vals_holder = zeros(length(species_counts), k_items);
bt_std_vals_holder = zeros(length(species_counts), k_items);
luke_std_vals_holder = zeros(length(species_counts), k_items);
btdist_mean_vals_holder = zeros(length(species_counts), k_items);
btdist_std_vals_holder = zeros(length(species_counts), k_items);

A = full_A(1, :);

for m = 1:length(species_counts)

    s = species_counts(m);
    curve = full_A(m + 1, :);
    figure(2); hold on;
    plot(curve, 'LineWidth', 2);
    A = cat(1, A, curve);
    A_norm = normalize_columns(A);

    num_elems = k_items - s + 1;
    bt_search_val_holder = zeros(num_repeat, num_elems);
    btdist_search_val_holder = zeros(num_repeat, num_elems);
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
            [min_inv_indices, min_inv_val] = random_search(A,s + k - 1, num_iters);
            bt_search_val_holder(r, k) = min_inv_val;

            %[bt_combo, bt_norm] = BT_dist_algo_v3(A, s + k - 1, 300, delta,4,10000);
            %[bt_combo, bt_norm] = og_dist_v1(A, s + k - 1, 10000, 5,10000);
            [bt_combo, bt_norm] = repeat_bt_opt(A, s+ k-1, 10000, 5, 10000, 5);

            btdist_search_val_holder(r, k) = bt_norm;

            current_iteration = current_iteration + 1;
            waitbar(current_iteration / total_iterations, wb, sprintf('k = %d, rep = %d...', k, r));
        end
    end

    close(wb);

    bt_mean_vals_holder(m, 1:num_elems) = mean(bt_search_val_holder, 1);
    bt_std_vals_holder(m, 1:num_elems) = std(bt_search_val_holder, 0, 1);
    btdist_mean_vals_holder(m, 1:num_elems) = mean(btdist_search_val_holder, 1);
    btdist_std_vals_holder(m, 1:num_elems) = std(btdist_search_val_holder, 0, 1);
    luke_mean_vals_holder(m, 1:num_elems) = luke_search_val_holder;
end

figure;
hold on;
for i = 1:length(species_counts)
    s = species_counts(i);
    num_elems = k_items - s + 1;
    errorbar(s:k_items, bt_mean_vals_holder(i, 1:num_elems), bt_std_vals_holder(i, 1:num_elems), 'LineWidth', 2, 'DisplayName', 'BT Fix', 'Color', colors(i, :));
    plot(s:k_items, luke_mean_vals_holder(i, 1:num_elems), '--', 'DisplayName', 'Luke', 'LineWidth', 2, 'Color', colors(i, :));
    errorbar(s:k_items, btdist_mean_vals_holder(i, 1:num_elems), btdist_std_vals_holder(i, 1:num_elems), ':', 'LineWidth', 2, 'DisplayName', 'BT Dist', 'Color', colors(i, :));
end
legend;

save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
idx = randi(9999);
file_name = ['multi_plot_data_with_btdist_', num2str(idx), '.mat'];
save(fullfile(save_folder, file_name), 'species_counts', 'k_items', 'wavelengths', 'bt_mean_vals_holder', 'bt_std_vals_holder', 'luke_mean_vals_holder', 'btdist_mean_vals_holder', 'btdist_std_vals_holder', 'A');
disp(['Data saved: ', file_name]);
end