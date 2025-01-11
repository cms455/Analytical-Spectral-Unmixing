% Generate random polynomial matrix
num_rows = 2;   % Number of rows (molecules)
num_cols = 100;  % Number of columns (wavelengths)
Nt = 2;         % Number of columns to select
max_iter = 100;  % Maximum iterations for refinement

% Generate a random polynomial matrix for testing
[A, x_values] = generate_polynomial_matrix(num_rows, num_cols);
%{
load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);
num_rows = size(A,1);
num_cols = size(A,2);
Nt = num_rows;
max_iter = 100;
%}
start_val = 125 + shift;
% Apply Bourgain–Tzafriri algorithm
[selected_indices, submatrix] = bourgain_tzafriri(A, Nt, max_iter);

% Display results
disp('Matrix A (rows = polynomial curves):');
disp(A);

disp('Selected column indices:');
disp(selected_indices);

disp('Selected submatrix:');
disp(submatrix);

disp('Condition Number:')
disp(cond(submatrix))

% Plot polynomial curves and selected wavelengths
figure;
plot(x_values, A', '-o', 'DisplayName', 'Polynomial Curves');
hold on;
scatter(x_values(selected_indices), submatrix(1, :), 80, 'r', 'filled', 'DisplayName', 'Selected Columns');
legend('show');
xlabel('Wavelength');
ylabel('Absorption Coefficient');
title('Random Polynomial Curves and Selected Columns (Bourgain–Tzafriri)');
grid on;


