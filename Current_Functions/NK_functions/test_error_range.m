min_w = 680;
max_w = 970;
species_bool = [1, 1, 0, 0, 0];
num_points = 400;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);

num_cols = size(A, 2);

%% Run random sampling of Frobenius norms
fro_norms = zeros(1, num_iters);

for i = 1:num_iters
    indices = randi(num_cols, 1, k);
    while numel(unique(indices)) < k
        indices = randi(num_cols, 1, k);
    end
    submatrix = A(:, indices);
    fro_norms(i) = norm(pinv(submatrix), 'fro');
end

%% Plot histogram
figure;
histogram(fro_norms, 50);
xlabel('Frobenius Norm of Pseudoinverse');
ylabel('Frequency');
title(sprintf('Distribution of Frobenius Norms for k = %d (n = %d samples)', k, num_iters));
set(gca, 'FontSize', 14);cl

%% Plot range as sorted line
figure;
plot(sort(fro_norms), 'k', 'LineWidth', 1.5);
xlabel('Sorted Sample Index');
ylabel('Frobenius Norm of Pseudoinverse');
title(sprintf('Range of Frobenius Norms for k = %d', k));
set(gca, 'FontSize', 14);
