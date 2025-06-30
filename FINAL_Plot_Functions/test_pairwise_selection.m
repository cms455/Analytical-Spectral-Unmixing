
A = randn(100, 50);  % example matrix


min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 0, 0];
num_points = 120;

k = 5;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
%A = randn(5, 120);  % example matrix
pair_map =visualize_pairwise_selection(A, k, 200000);

[best_combo,best_norm]=repeat_bt_opt(A, k, 10000, 5, 10000, 5);

disp('Regular Distribution');
disp(best_combo);
disp('\n')
sprintf('%.8f', best_norm)

tic;
[best_genotype, best_fitness] = nk_column_selector(A, k, 1000);
time = toc;
disp('NK');
sprintf('%.8f', best_fitness)
disp(toc);