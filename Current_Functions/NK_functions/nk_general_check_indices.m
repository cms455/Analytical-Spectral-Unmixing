min_w = 680;
max_w = 970;
species_bool = [1, 1, 1, 1, 1];
num_points = 120;
wavelengths = linspace(min_w, max_w, num_points);
num_species = sum(species_bool);
%load_A = load('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/combined_spectra_v2.mat','combined_spectra');
%load_A = load_A.combined_spectra;
%full_A = load_A(:,2:end)';
A = build_absorption_matrix(min_w, max_w, species_bool, num_points);
species_counts = [2, 3, 4, 5];
k = 2;
num_nk_iters = 1000;
num_iters = 1000000;

[l_submatrix, l_indices] = luke_algorithm(A', k);

[nk_indices_binary, nk_val] = nk_column_selector(A, k, num_nk_iters);
nk_indices = find(nk_indices_binary==1);
%[min_inv_indices, min_inv_val] = random_search(A, k, num_iters);

disp(' Luke Indices: ')
disp(l_indices)
disp('\n NK Indices:')
disp(nk_indices)
%disp('\n Random Indices:')
%disp(min_inv_indices)