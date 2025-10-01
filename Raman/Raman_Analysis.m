% Load the reference spectra
load_data = load('/Users/calvinsmith/Bouma_lab/Raman_Project/Raman-in-the-eye/CS_notebooks/refs.mat');
A = load_data.refs;        % size: [n_wavenumbers x n_refs]
wn = load_data.wn_grid;    % size: [n_wavenumbers x 1]

k = 30;
num_iters = 1000;

[best_wn, best_fitness] = nk_column_selector(A', k, num_iters);
disp(best_wn);
disp(abs(best_fitness));

% --- Get the selected wavenumbers ---
selected_idx = find(best_wn == 1);    % indices of chosen wn
selected_vals = wn(selected_idx);     % actual wavenumber values

disp("Selected_Wavelengths:")
disp(selected_vals)

for i = 1:length(selected_vals)
    disp()
end
% --- Plot spectra ---
figure;
plot(wn, A, 'LineWidth', 1.5);
xlabel('Wavenumber (cm^{-1})');
ylabel('Normalized Intensity');
title('Reference Spectra with Selected Wavenumbers');
grid on;
hold on;

% --- Add vertical lines at chosen wavenumbers ---
for i = 1:length(selected_vals)
    xline(selected_vals(i), 'r--', 'LineWidth', 1.5, ...
          'Label', sprintf('%d', selected_vals(i)), ...
          'LabelVerticalAlignment', 'bottom', ...
          'LabelHorizontalAlignment', 'right');
end

legend({'Oleic','Saturated','Protein'}, 'Location','best');
hold off;
