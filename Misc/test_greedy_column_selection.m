% Set parameters
num_rows = 2;   % Number of rows (molecules)
num_cols = 100;  % Number of columns (wavelengths)
Nt = 2;         % Number of columns to select

% Generate random polynomial matrix
[A, x_values] = generate_polynomial_matrix(num_rows, num_cols);

% Perform greedy column selection
[selected_indices, submatrix] = greedy_column_selection(A, Nt);

% Display results
disp('Matrix A (rows = polynomial curves):');
disp(A);

disp('Selected column indices:');
disp(selected_indices);

disp('Selected submatrix:');
disp(submatrix);

% Plot the polynomial curves and selected wavelengths
figure;
plot(x_values, A', '-o', 'DisplayName', 'Polynomial Curves');
hold on;
scatter(x_values(selected_indices), submatrix(1, :), 80, 'r', 'filled', 'DisplayName', 'Selected Columns');
legend('show');
xlabel('Wavelength');
ylabel('Absorption Coefficient');
title('Random Polynomial Curves and Selected Columns');
grid on;

function [A, x_values] = generate_polynomial_matrix(num_rows, num_cols)
    % Input:
    % num_rows: Number of rows (e.g., molecules)
    % num_cols: Number of columns (e.g., wavelengths)
    %
    % Output:
    % A: Matrix where each row represents a degree-3 polynomial curve
    % x_values: x-axis points (wavelengths) where the polynomials are sampled

    % Generate x-axis values (e.g., wavelengths)
    x_values = linspace(0, 5, num_cols);  % Sampled from 0 to 1

    % Initialize matrix
    A = zeros(num_rows, num_cols);

    % Generate random degree-3 polynomial coefficients for each row
    for i = 1:num_rows
        % Random coefficients for degree-3 polynomial: ax^3 + bx^2 + cx + d
        coeffs = rand(1, 4) - 0.5;  % Coefficients in the range [-0.5, 0.5]

        % Evaluate polynomial at x_values
        A(i, :) = polyval(coeffs, x_values) + 5;
    end
end