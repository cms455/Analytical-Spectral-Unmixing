min_w = 680;
max_w = 800;
species_bool = [1, 1, 1,1,1];
num_points = 60;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
full_A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
species_counts = [2,3, 4,5];

k_items = 8;

num_iters = 10000;
figure(1);
figure(2);
hold on;
for i = 1:min(species_counts)-1
    plot(full_A(i,:),'LineWidth',2);
    A = full_A(i,:);
end


bt_mean_vals_holder = zeros(length(species_counts),k_items);
luke_mean_vals_holder = zeros(length(species_counts),k_items);

bt_std_vals_holder = zeros(length(species_counts),k_items);
luke_std_vals_holder = zeros(length(species_counts),k_items);

for m = 1:length(species_counts)
s = species_counts(m);
curve = full_A(m+1,:);
figure(2);
hold on;
plot(curve, 'LineWidth',2);
A = cat(1,A,curve);
A_norm = normalize_columns(A);

num_elems = k_items -s + 1;

% Storage for results
bt_search_val_holder = zeros(num_repeat, num_elems);
bt_times = zeros(num_repeat, num_elems);

% Luke's algorithm storage
luke_search_val_holder = zeros(1, num_elems);
luke_times = zeros(1, num_elems);

% Initialize waitbar
total_iterations = num_repeat * (num_elems); % Total iterations for Bourgain-Tzafriri
current_iteration = 0; % Counter for current progress
wb = waitbar(0, 'Running algorithms, please wait...');

% Precompute Luke's algorithm results for all k
for k = 1:num_elems
    % Run Luke's algorithm once per k
    [l_submatrix, ~] = luke_algorithm(A', s+k-1);
    luke_search_val_holder(k) = norm(pinv(l_submatrix'), 'Fro');
end

% Run Bourgain-Tzafriri algorithm and update waitbar
for k = 1:num_elems
    for r = 1:num_repeat
        % Bourgain-Tzafriri Algorithm
        [~, min_inv_indices, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, s+k-1, num_iters);
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

bt_mean_vals_holder(m,1:num_elems) = bt_mean_vals;
bt_std_vals_holder(m,1:num_elems) = bt_std_vals;


luke_mean_vals = luke_search_val_holder; % No random component, so directly use results


luke_mean_vals_holder(m,1:num_elems) = luke_mean_vals;



end


figure;
hold on;
for i = 1:length(species_counts)
    s = species_counts(i);
    num_elems = k_items -s + 1;
    errorbar(s:k_items, bt_mean_vals_holder(i,1:num_elems), bt_std_vals_holder(i,1:num_elems),'LineWidth',2,'DisplayName','BT-ALGO', 'Color', colors(i, :));
    plot(s:k_items,luke_mean_vals_holder(i, 1:num_elems),'--','DisplayName','Luke Algorithm', 'LineWidth',2, 'Color', colors(i, :));
end
legend;




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

