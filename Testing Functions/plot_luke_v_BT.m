% Parameters
N = 100; % Number of iterations

% Containers for Bourgain-Tzafriri algorithm
condition_numbers_norm_BT = zeros(1, N); % Normalized condition numbers
times_BT = zeros(1, N); % Execution times

% Run Bourgain-Tzafriri Algorithm N times
for i = 1:N
    A = randn(10,1000);
    A_norm = normalize_columns(A);
    tic;
    [selected_indices, ~] = bourgain_tzafriri_v2(A_norm); % Run the Bourgain-Tzafriri algorithm
    times_BT(i) = toc;
    
    % Calculate condition number for normalized submatrix
    A_norm_selected = A_norm(:, selected_indices);
    condition_numbers_norm_BT(i) = cond(A_norm_selected);
end

% Plot normalized condition numbers for Bourgain-Tzafriri
mean_cond_norm_BT = mean(condition_numbers_norm_BT);
std_cond_norm_BT = std(condition_numbers_norm_BT);

figure;
histogram(condition_numbers_norm_BT, 20, 'FaceColor', [0, 0, 1], 'EdgeColor', 'black'); % Blue for normalized
hold on;
xline(mean_cond_norm_BT, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', mean_cond_norm_BT));
xline(mean_cond_norm_BT + std_cond_norm_BT, 'g--', 'LineWidth', 2, 'Label', sprintf('+1 Std Dev = %.2f', mean_cond_norm_BT + std_cond_norm_BT));
xline(mean_cond_norm_BT - std_cond_norm_BT, 'g--', 'LineWidth', 2, 'Label', sprintf('-1 Std Dev = %.2f', mean_cond_norm_BT - std_cond_norm_BT));
title('Bourgain-Tzafriri Algorithm: Normalized Condition Numbers');
xlabel('Condition Number');
ylabel('Frequency');
legend('Condition Numbers', 'Mean', '+/- 1 Std Dev');
hold off;

% Plot execution times for Bourgain-Tzafriri
mean_time_BT = mean(times_BT);
std_time_BT = std(times_BT);

figure;
histogram(times_BT, 20, 'FaceColor', [0, 0, 1], 'EdgeColor', 'black'); % Blue for times
hold on;
xline(mean_time_BT, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.4f s', mean_time_BT));
xline(mean_time_BT + std_time_BT, 'g--', 'LineWidth', 2, 'Label', sprintf('+1 Std Dev = %.4f s', mean_time_BT + std_time_BT));
xline(mean_time_BT - std_time_BT, 'g--', 'LineWidth', 2, 'Label', sprintf('-1 Std Dev = %.4f s', mean_time_BT - std_time_BT));
title('Bourgain-Tzafriri Algorithm: Execution Times');
xlabel('Time (s)');
ylabel('Frequency');
legend('Execution Times', 'Mean', '+/- 1 Std Dev');
hold off;

% Containers for Luke Algorithm
condition_numbers_norm_Luke = zeros(1, N); % Normalized condition numbers
times_Luke = zeros(1, N); % Execution times

% Run Luke Algorithm N times
for i = 1:N
    A = randn(10,1000);
    A_norm = normalize_columns(A);
    tic;
    [l_submatrix, l_indices] = luke_algorithm(A', 3); % Run Luke Algorithm
    times_Luke(i) = toc;

    % Calculate condition number for normalized submatrix
    A_norm_selected_Luke = A_norm(:, l_indices);
    condition_numbers_norm_Luke(i) = cond(A_norm_selected_Luke);
end

% Plot normalized condition numbers for Luke Algorithm
mean_cond_norm_Luke = mean(condition_numbers_norm_Luke);
std_cond_norm_Luke = std(condition_numbers_norm_Luke);

figure;
histogram(condition_numbers_norm_Luke, 20, 'FaceColor', [1, 0.5, 0], 'EdgeColor', 'black'); % Orange for normalized
hold on;
xline(mean_cond_norm_Luke, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.2f', mean_cond_norm_Luke));
xline(mean_cond_norm_Luke + std_cond_norm_Luke, 'g--', 'LineWidth', 2, 'Label', sprintf('+1 Std Dev = %.2f', mean_cond_norm_Luke + std_cond_norm_Luke));
xline(mean_cond_norm_Luke - std_cond_norm_Luke, 'g--', 'LineWidth', 2, 'Label', sprintf('-1 Std Dev = %.2f', mean_cond_norm_Luke - std_cond_norm_Luke));
title('Luke Algorithm: Normalized Condition Numbers');
xlabel('Condition Number');
ylabel('Frequency');
legend('Condition Numbers', 'Mean', '+/- 1 Std Dev');
hold off;

% Plot execution times for Luke Algorithm
mean_time_Luke = mean(times_Luke);
std_time_Luke = std(times_Luke);

figure;
histogram(times_Luke, 20, 'FaceColor', [1, 0.5, 0], 'EdgeColor', 'black'); % Orange for times
hold on;
xline(mean_time_Luke, 'r--', 'LineWidth', 2, 'Label', sprintf('Mean = %.4f s', mean_time_Luke));
xline(mean_time_Luke + std_time_Luke, 'g--', 'LineWidth', 2, 'Label', sprintf('+1 Std Dev = %.4f s', mean_time_Luke + std_time_Luke));
xline(mean_time_Luke - std_time_Luke, 'g--', 'LineWidth', 2, 'Label', sprintf('-1 Std Dev = %.4f s', mean_time_Luke - std_time_Luke));
title('Luke Algorithm: Execution Times');
xlabel('Time (s)');
ylabel('Frequency');
legend('Execution Times', 'Mean', '+/- 1 Std Dev');
hold off;
