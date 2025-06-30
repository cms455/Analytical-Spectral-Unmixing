
min_w = 680;
max_w = 970;
species_bool = [1, 1, 0, 0, 0];
num_points = 290;
%deltas = [4,6,10,14];
deltas = [100,100,200,200];
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
%full_A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
%load_A = load('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/combined_spectra_v2.mat','combined_spectra');
%load_A = load_A.combined_spectra;
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
%A = load_A(:,2:end)';
%A = A(1:4,:);
species_counts = [2, 3, 4, 5];
k = 2;
num_repeat = 1;
num_iters = 10000;


[l_submatrix, l_indices] = luke_algorithm(A', k);
luke_norm = norm(pinv(l_submatrix'), 'fro');
fprintf('\nLuke Norm: %d \n',luke_norm);
disp("Indices:")
disp(l_indices);

[combo, norm_val] = repeat_bt_opt(A, k, 20000, 5, 10000, 5);
fprintf('\nBest Norm: %.8f\n', norm_val);
disp('Best Combo:'); disp(combo);


%{
[rand_combo, rand_norm] = random_search(A, k, 1000000);
fprintf('\nRandom Norm: %d',rand_norm);
disp("\n Indices:")
disp(rand_combo);
%}