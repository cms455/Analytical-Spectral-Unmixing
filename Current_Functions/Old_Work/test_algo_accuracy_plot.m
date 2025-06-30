% Load Hb and HbO spectrum matrix
load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR_Spectrum.csv');
load_A = load_A()';
load_A = load_A(:,225:end);
% Extract actual wavelengths and spectral data
wavelengths = load_A(1, :); % Extract wavelength values from the second row
A = load_A(2:3,:)'; % Only Hb and HbO (assuming they are in rows 3 and 4)
A = A';
% Trim the data to remove unwanted shifts
A_norm = normalize_columns(A);

% Parameters
num_iters = 5000; % Number of iterations for BT algorithm
pick_cols = 2; % Number of selected wavelengths

% ---- Run Bourgain-Tzafriri Algorithm ----
disp('Running BT Algorithm...');
[~, bt_selected_indices, ~, ~, ~, ~] = bourgain_tzafriri_all_fix_selections(A, A_norm, pick_cols, num_iters);
bt_selected_wavelengths = wavelengths(bt_selected_indices);

% ---- Run Luke's Algorithm ----
disp('Running Luke Algorithm...');
[~, luke_selected_indices] = luke_algorithm(A', pick_cols);
luke_selected_wavelengths = wavelengths(luke_selected_indices);

% ---- Plot Spectral Curves with Selected Wavelengths ----
figure;
hold on;

% Plot Hb and HbO spectral curves
plot(wavelengths, A(1, :), 'r', 'LineWidth', 2, 'DisplayName', 'Hb');
plot(wavelengths, A(2, :), 'b', 'LineWidth', 2, 'DisplayName', 'HbO');

% Overlay BT-selected wavelengths as red dashed vertical lines
for wl = bt_selected_wavelengths
    xline(wl, 'r--', 'LineWidth', 1.5, 'Label', 'BT Selection');
end

% Overlay Luke-selected wavelengths as blue solid vertical lines
for wl = luke_selected_wavelengths
    xline(wl, 'b-', 'LineWidth', 1.5, 'Label', 'Luke Selection');
end

% Formatting
xlabel("Wavelength (nm)");
ylabel("Absorption");
title("Hb vs. HbO Spectral Curves with Selected Wavelengths");
legend('Location', 'Best');
grid on;
hold off;

% ---- Save Data for Future Use ----
save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
file_name = fullfile(save_folder, 'BT_Luke_Hb_HbO_selections.mat');

save(file_name, 'bt_selected_wavelengths', 'luke_selected_wavelengths', 'wavelengths', 'A');
disp(['Selection data saved as: ', file_name]);
