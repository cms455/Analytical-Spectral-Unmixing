% Load matrix
load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR.csv');
A = load_A(3:end, 2:end)';
shift = 225;
A = A(:, shift:end);
A_norm = normalize_columns(A);

% Parameters
max_runs = 2000;
num_iters = 100;
num_inner_iters = 5;
pick_cols = 2;
start_runs = 30;

% Initialize holders
time_holder = zeros(num_iters, num_inner_iters);
inv_val_holder = zeros(num_iters, num_inner_iters);
runs_list = linspace(start_runs, max_runs, num_iters);

 true_min_inv = 0.001137095318943;

% Initialize waitbar
total_steps = num_iters * num_inner_iters;
step = 0;
h = waitbar(0, 'Processing iterations...');

for k = 1:num_iters
    runs = round(runs_list(k));
    for j = 1:num_inner_iters
        % Simulate matrix and normalize
        %A = randn(10, 1000);
        %A_norm = normalize_columns(A);  % Assume normalize_columns is defined
        
        % Timer
        tic;
        [~, ~, ~, ~, ~, min_inv_val] = bourgain_tzafriri_k_cols_n_iters_EMD(A, A_norm, pick_cols, runs); % Function call
        time = toc;
        
        % Store values
        inv_val_holder(k, j) = min_inv_val;
        time_holder(k, j) = time;

        % Update waitbar
        step = step + 1;
        waitbar(step / total_steps, h, sprintf('Processing step %d of %d...', step, total_steps));
    end
end

% Close waitbar
close(h);



tic;
[l_submatrix,l_indices] = luke_algorithm(A',pick_cols);
toc;

luke_min_inv = norm(pinv(l_submatrix),'Fro');

% Calculate mean and standard deviation
time_mean = mean(time_holder, 2);
time_std = std(time_holder, 0, 2);
inv_val_mean = mean(inv_val_holder, 2);
inv_val_std = std(inv_val_holder, 0, 2);

% Plot time
figure;
errorbar(runs_list, time_mean, time_std, '-o');
xlabel('Number of Iterations');
ylabel('Time (s)');
title('Mean and Standard Deviation of Time');
grid on;

% Plot minimum inverse
figure;
errorbar(runs_list, inv_val_mean, inv_val_std, '-o');
yline(luke_min_inv, 'LineStyle',  '-','Color', 'r','Label','Luke Min Inverse');
yline(true_min_inv,'LineStyle',  '-','Color', 'g','Label','Luke Min Inverse')
xlabel('Number of Iterations');
ylabel('Norm of Inverse');
title('Mean and Standard Deviation of Condition Number');
grid on;

% Plot minimum inverse number with log-log scale
figure;
errorbar(runs_list, inv_val_mean, inv_val_std, '-o');
hold on;
yline(luke_min_inv, 'LineStyle', '-', 'Color', 'r', 'Label', 'Luke Min Inverse');
yline(true_min_inv, 'LineStyle', '-', 'Color', 'g', 'Label', 'True Min Inverse');
set(gca, 'YScale', 'log'); % Set log-log scale
xlabel('Number of Iterations (log scale)');
ylabel('Norm of Inverse (log scale)');
%ylim([true_min_inv -(1e-6),10e-1])
title('Mean and Standard Deviation of Minimum Inverse (Log-Log)');
grid on;
legend('Mean ± STD', 'Luke Min Inverse', 'True Min Inverse');
hold off;
