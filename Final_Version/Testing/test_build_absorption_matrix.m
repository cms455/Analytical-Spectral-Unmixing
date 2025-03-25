min_w = 600;
max_w = 1000;
species_bool = [1,1,1,0];
num_points = 40;
wavelengths = linspace(min_w,max_w, num_points);

A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
figure;
hold on;
for i =  1:size(A,1)
    plot(A(i,:));
end
legend;

A_norm = normalize_columns(A);
k_items = 5;
num_iters = 50000;
num_species = sum(species_bool);

% Storage for results
bt_search_val_holder = zeros(1, k_items - 1);
luke_search_val_holder = zeros(1, k_items - 1);
bt_times = zeros(1, k_items - 1);
luke_times = zeros(1, k_items - 1);

for k = num_species:k_items
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

    % Plot all rows of A with selected indices
    figure;
    hold on;
    for i = 1:num_species
        plot(wavelengths, A(i, :), 'LineWidth', 2); % Plot all rows
    end

    % Add red lines for BT selected indices and blue lines for Luke selected indices
    for i = 1:k
        x_value_bt = wavelengths(min_inv_indices(i));
        xline(x_value_bt, 'r--', 'LineWidth', 1.5); % Red dashed lines for BT
        x_value_luke = wavelengths(l_indices(i));
        xline(x_value_luke, 'b--', 'LineWidth', 1.5); % Blue dashed lines for Luke
    end
    hold off;

    % Customize the plot
    xlabel('Wavelength (nm)');
    ylabel('Absorption');
    title(sprintf('k = %d: Absorption Spectra with Selected Indices Highlighted', k));
    legend(arrayfun(@(x) sprintf('Species %d', x), 1:num_species, 'UniformOutput', false), 'Location', 'Best');
  
end

% Plot minimum inverse values comparison
figure;
hold on;
set(gca,'FontSize',14)
plot(num_species:k_items, bt_search_val_holder(num_species-1:end), '-o', 'DisplayName', 'Bourgain-Tzafriri', 'LineWidth',2);
plot(num_species:k_items, luke_search_val_holder(num_species-1:end), '-o', 'DisplayName', 'Luke Algorithm', 'LineWidth',2);
hold off;
xlabel('k (Wavelength Selections)');
ylabel('Minimum Inverse Frobenius Norm');
title('Comparison of Minimum Inverse Values for BT and Luke Algorithms');
legend('Location', 'Best');


% Plot runtime comparison
figure;
hold on;
set(gca,'FontSize',14)
plot(num_species:k_items, bt_times(num_species-1:end), '-o', 'DisplayName', 'Bourgain-Tzafriri', 'LineWidth',2);
plot(num_species:k_items, luke_times(num_species-1:end), '-o', 'DisplayName', 'Luke Algorithm', 'LineWidth',2);
hold off;
xlabel('k (Wavelength Selections)');
ylabel('Runtime (seconds)');
title('Comparison of Runtime for BT and Luke Algorithms');
legend('Location', 'Best');



