% Generate random polynomial matrix
num_rows = 2;   % Number of rows (molecules)
num_cols = 100;  % Number of columns (wavelengths)
Nt = 2;         % Number of columns to select
max_iter = 100;  % Maximum iterations for refinement

% Generate a random polynomial matrix for testing
%[A, x_values] = generate_polynomial_matrix(num_rows, num_cols);

load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end);
A = A';
shift = 225;
A = A(:, shift:end);

num_rows = size(A,1);
num_cols = size(A,2);
Nt = 4;
max_iter = 100;

%A = randn(10,1000);
%A = normalize_columns(A);
start_val = 125 + shift;
% Apply Bourgain–Tzafriri algorithm
tic;
[selected_indices, submatrix] = bourgain_tzafriri(A, Nt, max_iter);
toc;
% Display results
%disp('Matrix A (rows = polynomial curves):');
%disp(A);

%disp('Selected column indices:');
%disp(selected_indices);

%disp('Selected submatrix:');%
%disp(submatrix);

disp('Condition Number:')
disp(cond(submatrix))

x_values = 2*((start_val):(start_val+num_cols-1));
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


% Build the noise matrix
%{
N = randn(num_rows,1);

C_n = submatrix .* N;

mean_C_n = sum(C_n,'all')/(num_rows*Nt);

std_C_n = std(C_n);

max_C_n = max(C_n,[],'all');

disp('C_n:')
disp(round(C_n))
fprintf('Mean C_n: %d\n', mean_C_n);
%}

idx = [1,144];
test = A(:,idx);
test_inv = pinv(test);

idx_2 = [1,135];
test_2 = A(:,idx_2);
test_2_inv = pinv(test_2);

submatrix_inv = pinv(submatrix);

fprintf('Condition Number of TEST: %d \n', cond(test));
fprintf('Condition Number of TEST2: %d \n', cond(test_2_inv));
fprintf('Selected Indices of TEST: %d \n', idx);
fprintf('Selected Indices of TEST 2: %d \n', idx_2);
%{
analyze_submatrix_with_noise(test_inv, 100);
analyze_submatrix_with_noise(submatrix_inv, 100);
%}

num_noise_distributions = 100;  % Number of noise distributions
[test_row_means, test_row_stds, test_total_mean, test_total_std] = analyze_submatrix_with_noise(test_inv, num_noise_distributions);
[submatrix_row_means, submatrix_row_stds, submatrix_total_mean, submatrix_total_std] = analyze_submatrix_with_noise(test_2_inv, num_noise_distributions);

% Plot comparison of row-wise statistics
figure;

% Row-wise means
subplot(2, 2, 1);
plot(1:length(test_row_means), test_row_means, '-o', 'DisplayName', 'test\_inv');
hold on;
plot(1:length(submatrix_row_means), submatrix_row_means, '-x', 'DisplayName', 'submatrix\_inv');
title('Row-wise Means');
xlabel('Row Index');
ylabel('Mean');
legend('show');
grid on;

% Row-wise standard deviations
subplot(2, 2, 2);
plot(1:length(test_row_stds), test_row_stds, '-o', 'DisplayName', 'test\_inv');
hold on;
plot(1:length(submatrix_row_stds), submatrix_row_stds, '-x', 'DisplayName', 'submatrix\_inv');
title('Row-wise Standard Deviations');
xlabel('Row Index');
ylabel('Standard Deviation');
legend('show');
grid on;

% Total mean comparison
subplot(2, 2, 3);
bar([test_total_mean, submatrix_total_mean]);
set(gca, 'XTickLabel', {'test\_inv', 'submatrix\_inv'});
title('Total Mean');
ylabel('Mean');
grid on;

% Total standard deviation comparison
subplot(2, 2, 4);
bar([test_total_std, submatrix_total_std]);
set(gca, 'XTickLabel', {'test\_inv', 'submatrix\_inv'});
title('Total Standard Deviation');
ylabel('Standard Deviation');
grid on;

% Enhance layout
sgtitle('Comparison of Noise Analysis Results');
