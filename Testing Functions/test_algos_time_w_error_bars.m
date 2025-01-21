% Load and preprocess the matrix A
load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end)';
shift = 225;
A = A(:, shift:end);

% Parameters
num_rows = size(A, 1);
num_cols = size(A, 2);

max_cols = 10;
num_iters = 10;
j = 3; % Number of repetitions for each column count
time_holder = zeros(j, num_iters);
cols = linspace(2, max_cols, num_iters);





% Main loop: Iterate over different numbers of columns
for k = 1:num_iters
    k_col = round(cols(k));
    A_rand = randn(10,200); % Generate random matrix for testing
    A_norm = normalize_columns(A_rand); % Normalize columns
    for rep = 1:j
        tic;
        bourgain_tzafriri_all_fix_selections(A_rand, A_norm, k_col, 500);
        time_holder(rep, k) = toc; % Store elapsed time for this run
    end
end

% Compute mean and standard deviation of times for error bars
mean_times = mean(time_holder, 1);
std_times = std(time_holder, 0, 1);

% Plot the results with error bars
figure;
errorbar(cols, mean_times, std_times, 'LineWidth', 2, 'CapSize', 8);
title('Time for Mod-BT Algorithm with Error Bars');
%ylim([0,0.05])
xlabel('Number of Columns');
ylabel('Time (s)');
set(gca, 'FontSize', 14);
grid on;
