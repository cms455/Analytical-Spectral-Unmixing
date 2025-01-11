% Parameters
num_iterations = 1000;
matrix_size = 10; % E is a square matrix
noise_level = 1;

% Initialize arrays for storing metrics
frobenius_norms = zeros(1, num_iterations);
least_singular_values = zeros(1, num_iterations);
C_n_values = zeros(1, num_iterations);

% Perform 100 iterations
for i = 1:num_iterations
    % Generate random matrix E and noise vector N
    E = randn(2,matrix_size);
    N = noise_level * randn(matrix_size, 1);
    
    % Calculate C_n = sqrt(sum(C_n.^2)), where C_n = E * N
    C_n = E * N;
    C_n_values(i) = sqrt(sum(C_n.^2));
    
    % Calculate Frobenius norm and least singular value of E
    frobenius_norms(i) = norm(E, 'fro');
    singular_values = svd(E);
    least_singular_values(i) = singular_values(end);
end

% Plot 1: C_n vs Least Singular Value
figure;
scatter(least_singular_values, C_n_values, 'filled');
xlabel('Least Singular Value of E');
ylabel('C_n Value');
title('C_n vs Least Singular Value of E');
grid on;

% Plot 2: C_n vs Frobenius Norm
figure;
scatter(frobenius_norms, C_n_values, 'filled');
xlabel('Frobenius Norm of E');
ylabel('C_n Value');
title('C_n vs Frobenius Norm of E');
grid on;
