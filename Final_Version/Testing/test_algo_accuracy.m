% Load matrix
load_A = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/HbOHbR_Spectrum.csv');
A = load_A(3:end, 2:end)';
shift = 225;
A = A(:, shift:end);
A_norm = normalize_columns(A);

% Parameters
max_runs = 5000;
num_iters = 20;
num_inner_iters = 4;
pick_cols = 2 ;
start_runs = 30;

% Initialize holders
time_holder = zeros(num_iters, num_inner_iters);
inv_val_holder = zeros(num_iters, num_inner_iters);
runs_list = linspace(start_runs, max_runs, num_iters);

% true_min_inv = 9.2557e-4;
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
        %[~, ~, ~, ~, ~, min_inv_val] = bourgain_tzafriri_k_cols_n_iters_EMD(A, A_norm, pick_cols, runs); % Function call
        [~, ~, ~, ~, ~, min_inv_val] = bourgain_tzafriri_all_fix_selections(A, A_norm, pick_cols, runs); % Function call
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
errorbar(runs_list, time_mean, time_std, '-o', 'LineWidth',2);
xlabel('Number of Iterations');
ylabel('Time (s)');
%title('Time for Modified BT');
grid on;
set(gca,'FontSize',14);


% Plot minimum inverse number with log-log scale
runs_list(runs_list < 0) = 0;
figure;
ax = axes;
plot(1:num_iters,inv_val_mean);
ylim([0,0.00124]);
ax.YScale = 'log';
drawnow;
ax.YLim  



figure;
ax = axes;
errorbar(1:num_iters, inv_val_mean, inv_val_std, '-o', 'LineWidth',2);
hold on;
yline(luke_min_inv,  'LineWidth',2,'LineStyle', '-', 'Color', 'r', 'Label', 'Greedy Algorithm');
yline(true_min_inv, 'LineWidth',2,'LineStyle', '-', 'Color', 'g', 'Label', 'Optimal Norm');
set(gca, 'FontSize',14);
set(gca, 'YScale', 'log');
  % Customize as per your data
ax.YScale = 'log';
ylim([0, max(inv_val_mean) * 1.1]);
xlabel('Number of Iterations');
ylabel('Norm of Inverse');
%title('Mean and Standard Deviation of Minimum Inverse ');
%legend('Mean of Mod-BT', 'Greedy Algo', 'Optimal Norm');
hold off;

% Save all necessary data for recreating the figures
save_folder = '/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/ASU_plot_data';

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

idx = randi(9999);

file_name = ['test_algo_accuracy_data_', num2str(idx),'_'];

% Save comparison data
save(fullfile(save_folder, file_name), ...
    'A', 'A_norm', 'shift', ...
    'max_runs', 'num_iters', 'num_inner_iters', 'pick_cols', 'start_runs', ...
    'runs_list', 'time_holder', 'inv_val_holder', ...
    'true_min_inv', 'luke_min_inv', ...
    'time_mean', 'time_std', 'inv_val_mean', 'inv_val_std');

disp('All necessary data for recreating the figures has been saved.');
