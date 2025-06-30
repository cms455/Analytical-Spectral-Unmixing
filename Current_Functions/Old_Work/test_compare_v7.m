min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 1, 1]; % Five species
num_points = 120;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);

% Plot spectral absorption curves
figure;
hold on;
for i = 1:num_species
    plot(wavelengths, A(i,:), 'LineWidth', 2);
end
xlabel('Wavelength (nm)');
ylabel('Absorption');
title('Spectral Absorption Curves');
grid on;

A_norm = normalize_columns(A);
k_items = 5;
num_iters = 10000; % Number of iterations for BT Algorithm

% Run Luke's Algorithm on the full curve
[l_submatrix, luke_selected_indices] = luke_algorithm(A', k_items);
luke_inverse_value = norm(pinv(l_submatrix), 'Fro');

% Run BT Algorithm once (no repetitions)
[~, bt_selected_indices, ~, ~, ~, bt_inverse_value] = bourgain_tzafriri_all_fix_selections(A, A_norm, k_items, num_iters);

% Overlay selected wavelengths on spectral curves
selected_luke_wavelengths = wavelengths(luke_selected_indices);
selected_bt_wavelengths = wavelengths(bt_selected_indices);

% Plot selected wavelengths
for w = 1:length(selected_luke_wavelengths)
    xline(selected_luke_wavelengths(w), 'b', 'LineWidth', 2, 'Label', 'Luke');
end
for w = 1:length(selected_bt_wavelengths)
    xline(selected_bt_wavelengths(w), 'r', 'LineWidth', 2, 'Label', 'BT');
end
legend('Location', 'Best');
hold off;

% Plot minimum inverse values for Luke and BT algorithms
figure;
hold on;
plot(k_items, bt_inverse_value, 'ro', 'MarkerSize', 10, 'DisplayName', 'BT Algorithm');
plot(k_items, luke_inverse_value, 'bo', 'MarkerSize', 10, 'DisplayName', 'Luke Algorithm');
xlabel('k (Wavelength Selections)');
ylabel('Minimum Inverse Frobenius Norm');
title('Comparison of Minimum Inverse Values for BT and Luke Algorithms');
legend('Location', 'Best');
grid on;
hold off;
