min_w = 680;
max_w = 970;
num_points = 290;
species_bool = [1, 1, 0, 0, 0];
num_species = sum(species_bool);

% Generate wavelengths
wavelengths = linspace(min_w, max_w, num_points);

% Build absorption matrix
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);

% Preallocate for range 40:end
start_idx = 40;
valid_range = start_idx:num_points;
num_eval = length(valid_range);

cond_vals = zeros(num_eval, 1);
pinv_norm_vals = zeros(num_eval, 1);

% Loop over all wavelengths to pair with the first column
for j = 1:num_eval
    i = valid_range(j);
    subA = A(:, [1, i]);  % Always column 1 + current index

    cond_vals(j) = cond(subA);
    pinv_norm_vals(j) = norm(pinv(subA), 'fro');
end

% Compute combined metric
combined_metric = (1/100)*(cond_vals) .* pinv_norm_vals;
eval_wavelengths = wavelengths(valid_range);

% Plot condition number
figure;
plot(eval_wavelengths, cond_vals, 'b-', 'LineWidth', 2);
xlabel('Second Wavelength (nm)');
ylabel('Condition Number');
title('Condition Number vs Second Wavelength');
grid on;

% Plot Frobenius norm
figure;
plot(eval_wavelengths, pinv_norm_vals, 'r--', 'LineWidth', 2);
xlabel('Second Wavelength (nm)');
ylabel('Frobenius Norm of Pseudoinverse');
title('Frobenius Norm of Pseudoinverse vs Second Wavelength');
grid on;

% Plot combined metric
figure;
plot(eval_wavelengths, combined_metric, 'm-', 'LineWidth', 2);
xlabel('Second Wavelength (nm)');
ylabel('Condition Number × Frobenius Norm');
title('Combined Metric vs Second Wavelength');
grid on;

% Find minima
[~, min_cond_idx] = min(cond_vals);
[~, min_pinv_idx] = min(pinv_norm_vals);
[~, min_combined_idx] = min(combined_metric);

min_cond_wavelength = eval_wavelengths(min_cond_idx);
min_pinv_wavelength = eval_wavelengths(min_pinv_idx);
min_combined_wavelength = eval_wavelengths(min_combined_idx);

% Print results with high precision
fprintf('Wavelength minimizing condition number: %.8f nm\n', min_cond_wavelength);
fprintf('Wavelength minimizing Frobenius norm of pseudoinverse: %.8f nm\n', min_pinv_wavelength);
fprintf('Wavelength minimizing combined metric: %.8f nm\n', min_combined_wavelength);
