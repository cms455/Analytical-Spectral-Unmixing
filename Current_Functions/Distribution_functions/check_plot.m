min_w = 680;
max_w = 970;
num_points = 290;
species_bool = [1, 1, 0, 0, 0];
num_species = sum(species_bool);

% Generate wavelengths
wavelengths = linspace(min_w, max_w, num_points);

% Build absorption matrix
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);

% Initialize results
num_windows = num_points - num_species + 1;
cond_vals = zeros(num_windows, 1);
pinv_norm_vals = zeros(num_windows, 1);
window_centers = zeros(num_windows, 1);

% Loop over sliding windows of consecutive wavelengths
for i = 1:num_windows
    idx = i:(i + num_species - 1);
    subA = A(:, idx);
    
    % Condition number and Frobenius norm of pseudoinverse
    cond_vals(i) = cond(subA);
    pinv_norm_vals(i) = norm(pinv(subA), 'fro');
    
    % Record center wavelength for plotting
    window_centers(i) = mean(wavelengths(idx));
end

% Plotting
figure;
yyaxis left;
plot(window_centers, cond_vals, 'b-', 'LineWidth', 2);
ylabel('Condition Number');
xlabel('Center Wavelength (nm)');
title('Condition Number and Pseudoinverse Frobenius Norm vs Wavelength Window');

yyaxis right;
plot(window_centers, pinv_norm_vals, 'r--', 'LineWidth', 2);
ylabel('Frobenius Norm of Pseudoinverse');

legend('Condition Number', 'Pseudoinverse Frobenius Norm');
grid on;
