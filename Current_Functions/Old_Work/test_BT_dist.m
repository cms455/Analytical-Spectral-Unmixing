min_w = 680;
max_w = 1000;
species_bool = [1, 1, 0, 0, 0];
num_points = 150;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);

delta = 0.8;
num_iters = 1000;
num_columns = size(A, 2);
selection_counts = zeros(1, num_columns);

max_time = 3; % seconds

h = waitbar(0, 'Running iterations...');

for i = 1:num_iters
    waitbar(i / num_iters, h, sprintf('Iteration %d of %d', i, num_iters));
    try
        tic;
        [conditioned_indices, submatrix] = bourgain_tzafriri_v2(A, delta);
        elapsed_time = toc;

        if elapsed_time > max_time
            fprintf("Iteration %d skipped (%.2f sec > %.2f sec)\n", i, elapsed_time, max_time);
            continue;
        end

        selection_counts(conditioned_indices) = selection_counts(conditioned_indices) + 1;

        % Optional: Uncomment for diagnostics
        % k_A = cond(A);
        % nK_A = cond(A(:,conditioned_indices));
        % nE_A = norm(pinv(A(:,conditioned_indices)),'Fro');
        % fprintf("\nIteration %d\n", i);
        % fprintf("OG K: %f \n", k_A);
        % fprintf("New K: %f \n", nK_A);
        % fprintf("New E: %f \n", nE_A);
        % disp("Indices:");
        % disp(conditioned_indices);

    catch ME
        fprintf("Iteration %d failed: %s\n", i, ME.message);
    end
end

close(h);

% Plot distribution
figure;
bar(selection_counts, 'FaceColor', [0.2 0.4 0.6]);
xlabel('Column Index');
ylabel('Selection Frequency');
title(sprintf('Column Selection Distribution over %d Iterations (\\delta = %.1f)', num_iters, delta));
grid on;


