
min_w = 680;
max_w = 800;
species_bool = [1, 1, 1, 1, 1];
num_points = 120;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
load_A = load('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/combined_spectra_v2.mat','combined_spectra');
load_A = load_A.combined_spectra;
%full_A = load_A(:,2:end)';
full_A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
A = full_A

%A = randn(50,100);
k = 10;
num_iters = 1000;
T_init = 1.0;
T_min = 1e-6;
alpha = 0.98;
kB = 1;

tic; 
[genotype, fitness] = nk_column_selector_anneal(A, k, num_iters, T_init, T_min, alpha, kB);
anneal_time = toc;
disp(find(genotype));
fprintf('Best fitness: %.6f\n', fitness);
fprintf('Anneal Time: %.6f\n', anneal_time);

tic;
[genotype, fitness] = nk_column_selector(A, k, num_iters);
reg_time = toc;
disp(find(genotype));
fprintf('Best fitness: %.6f\n', fitness);
fprintf('Regular Time: %.6f\n', reg_time);